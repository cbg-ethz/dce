run.all.models <- function(
  wt.graph, wt.X,
  mt.graph, mt.X,
  wt.graph.perturbed,
  beta.magnitude,
  methods = NULL,
  effect.type,
  adjustment.type,
  latent
) {
  # for null model
  negweight.range <- c(-beta.magnitude, 0)
  posweight.range <- c(0, beta.magnitude)

  # ground truth
  if (effect.type == "direct") {
      # direct effects:
      ground.truth <- list(dce = as(mt.graph,"matrix") - as(wt.graph,"matrix"))
  } else if (effect.type == "total") {
      # total effects
      ground.truth <- list(dce = trueEffects(mt.graph) - trueEffects(wt.graph))
  }
  # reduce data for correlations
  graph <- as(wt.graph, "matrix")
  wt.X.cor <- wt.X[, seq_len(ncol(graph))]
  mt.X.cor <- mt.X[, seq_len(ncol(graph))]

  # correlations
  time.tmp <- Sys.time()
  if (is.null(methods) || "cor" %in% methods) {
    res.cor <- list(dce = cor(mt.X.cor) - cor(wt.X.cor))
    res.cor$dce_pvalue <- pcor_perm(wt.X.cor, mt.X.cor, fun = cor)
  } else {
    res.cor <- ground.truth
    res.cor$dce_pvalue <- ground.truth$dce*0
  }
  res.cor$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
  res.cor$dce_pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
  time.cor <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  time.tmp <- Sys.time()
  if (is.null(methods) || "pcor" %in% methods) {
    res.pcor <- list(dce = pcor(mt.X.cor) - pcor(wt.X.cor))
    res.pcor$dce_pvalue <- pcor_perm(wt.X.cor, mt.X.cor, fun = pcor)
  } else {
    res.pcor <- ground.truth
    res.pcor$dce_pvalue <- ground.truth$dce*0
  }
  res.pcor$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
  res.pcor$dce_pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
  time.pcor <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  # differential causal effects
  if (link.method == "log") {
    solver.args <- list(method = "glm.fit", link = "log")
  } else {
    solver.args <- list(method = "glm.dce.nb.fit", link = "identity")
  }

  time.tmp <- Sys.time()
  if (is.null(methods) || "dce" %in% methods) {
    res.dce <- dce::dce_nb(
      wt.graph.perturbed, wt.X, mt.X,
      adjustment_type = adjustment.type,
      effect_type = effect.type,
      solver_args = solver.args,
      lib_size = TRUE
    )
  } else {
    res.dce <- ground.truth
    res.dce$dce_pvalue <- ground.truth$dce*0
    res.dce$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
    res.dce$dce_pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
  }
  time.dce <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  time.tmp <- Sys.time()
  if (is.null(methods) || "dce.log" %in% methods) {
    solver.args.log <- list(method = "glm.fit", family = gaussian)
    res.dce.log <- dce::dce(
      wt.graph.perturbed, log(wt.X+1), log(mt.X+1),
      solver = "glm2",
      adjustment_type = adjustment.type,
      effect_type = effect.type,
      solver_args = solver.args.log,
      lib_size = TRUE
    )
  } else {
    res.dce.log <- ground.truth
    res.dce.log$dce_pvalue <- ground.truth$dce*0
    res.dce.log$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
    res.dce.log$dce_pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
  }
  time.dce.log <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  time.tmp <- Sys.time()
  if (is.null(methods) || "dce.latent" %in% methods) {
    res.dce.latent <- dce::dce_nb(
      wt.graph.perturbed, wt.X, mt.X,
      adjustment_type = adjustment.type,
      effect_type = effect.type,
      solver_args = solver.args,
      lib_size = TRUE,
      latent = latent
    )
  } else {
    res.dce.latent <- ground.truth
    res.dce.latent$dce_pvalue <- ground.truth$dce*0
    res.dce.latent$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
    res.dce.latent$dce_pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
  }
  time.dce.latent <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  time.tmp <- Sys.time()
  if (is.null(methods) || "dce.lr" %in% methods) {
    res.dce.lr <- dce::dce_nb(
      wt.graph.perturbed, wt.X, mt.X,
      adjustment_type = adjustment.type,
      effect_type = effect.type,
      solver_args = solver.args, test = "lr",
      lib_size = TRUE
    )
  } else {
    res.dce.lr <- ground.truth
    res.dce.lr$dce_pvalue <- ground.truth$dce*0
    res.dce.lr$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
    res.dce.lr$dce_pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
  }
  time.dce.lr <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  time.tmp <- Sys.time()
  if (is.null(methods) || "dce.nolib" %in% methods) {
    res.dce.nolib <- dce::dce_nb(
      wt.graph.perturbed, wt.X, mt.X,
      adjustment_type = adjustment.type,
      effect_type = effect.type,
      solver_args = solver.args,
      lib_size = FALSE
    )
  } else {
    res.dce.nolib <- ground.truth
    res.dce.nolib$dce_pvalue <- ground.truth$dce*0
    res.dce.nolib$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
    res.dce.nolib$dce_pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
  }
  time.dce.nolib <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  time.tmp <- Sys.time()
  if (is.null(methods) || "dce.tpm" %in% methods) {
    compute.tpm <- function(counts, len.kb = 1) {
      # TODO: is len.kb == 1 reasonable in our case?
      rpk <- counts / len.kb
      pm.scale <- rowSums(rpk) / 1e6
      return(rpk / pm.scale)
    }

    res.dce.tpm <- dce::dce_nb(
      wt.graph.perturbed, compute.tpm(wt.X), compute.tpm(mt.X),
      adjustment_type = adjustment.type,
      effect_type = effect.type,
      solver_args = solver.args,
      lib_size = FALSE
    )
  } else {
    res.dce.tpm <- ground.truth
    res.dce.tpm$dce_pvalue <- ground.truth$dce*0
    res.dce.tpm$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
    res.dce.tpm$dce_pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
  }
  time.dce.tpm <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  time.tmp <- Sys.time()
  if (is.null(methods) || "dce.tpmlog" %in% methods) {
    compute.tpm <- function(counts, len.kb = 1) {
      # TODO: is len.kb == 1 reasonable in our case?
      rpk <- counts / len.kb
      pm.scale <- rowSums(rpk) / 1e6
      return(rpk / pm.scale)
    }

    solver.args.log <- list(method = "glm.fit", family = gaussian)
    res.dce.tpmlog <- dce::dce(
      wt.graph.perturbed, log(compute.tpm(wt.X)+1), log(compute.tpm(mt.X)+1),
      solver = "glm2",
      adjustment_type = adjustment.type,
      effect_type = effect.type,
      solver_args = solver.args.log,
      lib_size = FALSE
    )
  } else {
    res.dce.tpmlog <- ground.truth
    res.dce.tpmlog$dce_pvalue <- ground.truth$dce*0
    res.dce.tpmlog$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
    res.dce.tpmlog$dce_pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
  }
  time.dce.tpmlog <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  time.tmp <- Sys.time()
  if (is.null(methods) || "dce.lm" %in% methods) {
    res.dce.lm <- dce::dce(
      wt.graph.perturbed, wt.X, mt.X,
      solver = "lm",
      adjustment_type = adjustment.type,
      effect_type = effect.type,
      lib_size = TRUE
    )
  } else {
    res.dce.lm <- ground.truth
    res.dce.lm$dce_pvalue <- ground.truth$dce*0
    res.dce.lm$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
    res.dce.lm$dce_pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
  }
  time.dce.lm <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  # null models
  time.tmp <- Sys.time()
  if (is.null(methods) || "rand" %in% methods) {
    tmp <- as(wt.graph.perturbed, "matrix") * NA
    tmp[which(as(wt.graph.perturbed, "matrix") != 0)] = (
      runif(sum(as(wt.graph.perturbed, "matrix") != 0), negweight.range[1], posweight.range[2]) -
        runif(sum(as(wt.graph.perturbed, "matrix") != 0), negweight.range[1], posweight.range[2])
    )
    tmp.pvals <- as(wt.graph.perturbed, "matrix") * NA
    tmp.pvals[which(as(wt.graph.perturbed, "matrix") != 0)] <- runif(sum(as(wt.graph.perturbed, "matrix") != 0), 0, 1)
    res.rand <- list(dce = tmp, dce_pvalue = tmp.pvals)
  } else {
    res.rand <- ground.truth
    res.rand$dce_pvalue <- ground.truth$dce*0
    res.rand$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
    res.rand$dce_pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
  }
  time.rand <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  # other methods
  time.tmp <- Sys.time()
  if (is.null(methods) || "causaldag" %in% methods) {
    # run causaldag
    dname.tmp <- "tmp.causaldag/"

    unlink(dname.tmp, recursive = TRUE)
    dir.create(dname.tmp, recursive = TRUE)

    write.csv(wt.X[, nodes(wt.graph.perturbed)], file.path(dname.tmp, "X_wt.csv"), quote=FALSE)
    write.csv(mt.X[, nodes(wt.graph.perturbed)], file.path(dname.tmp, "X_mt.csv"), quote=FALSE)

    cmd <- glue::glue("python3 execute_causaldag.py {dname.tmp}/X_wt.csv {dname.tmp}/X_mt.csv {dname.tmp}/difference_matrix.csv")
    #print(cmd)
    system(cmd)

    if (file.exists(file.path(dname.tmp, "difference_matrix.csv"))) {
      df.causaldag <- read.csv(file.path(dname.tmp, "difference_matrix.csv"))

      rownames(df.causaldag) <- df.causaldag$X
      df.causaldag <- subset(df.causaldag, select=-X)
    } else {
      df.causaldag <- NULL
    }

    unlink(dname.tmp, recursive = TRUE)

    # extract results
    if (is.null(df.causaldag)) {
      # execution crashed
      res.causaldag <- NULL # TODO: make this better
      print("NOOOO")
    } else {
      # everything went fine
      res.causaldag <- list(
        dce = as.matrix(df.causaldag),
        dce_pvalue = as.matrix(df.causaldag)
      )

      res.causaldag$dce_pvalue[res.causaldag$dce_pvalue != 0] <- 0.01 # something significant...
      res.causaldag$dce_pvalue[res.causaldag$dce_pvalue == 0] <- 0.99

      res.causaldag$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
      res.causaldag$dce_pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
    }
  } else {
    res.causaldag <- ground.truth
    res.causaldag$dce_pvalue <- ground.truth$dce*0
    res.causaldag$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
    res.causaldag$dce_pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
  }
  time.causaldag <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  # gather results
  df.res <- data.frame(
    truth=as.vector(ground.truth$dce),
    cor=as.vector(res.cor$dce),
    pcor=as.vector(res.pcor$dce),
    dce=as.vector(res.dce$dce),
    dce.log=as.vector(res.dce.log$dce),
    dce.latent=as.vector(res.dce.latent$dce),
    dce.lr=as.vector(res.dce.lr$dce),
    dce.nolib=as.vector(res.dce.nolib$dce),
    dce.tpm=as.vector(res.dce.tpm$dce),
    dce.tpmlog=as.vector(res.dce.tpmlog$dce),
    dce.lm=as.vector(res.dce.lm$dce),
    rand=as.vector(res.rand$dce),
    causaldag=as.vector(res.causaldag$dce)
  )

  df.pvalues <- data.frame(
    truth=as.vector(ground.truth$dce),
    cor=as.vector(res.cor$dce_pvalue),
    pcor=as.vector(res.pcor$dce_pvalue),
    dce=as.vector(res.dce$dce_pvalue),
    dce.log=as.vector(res.dce.log$dce_pvalue),
    dce.latent=as.vector(res.dce.latent$dce_pvalue),
    dce.lr=as.vector(res.dce.lr$dce_pvalue),
    dce.nolib=as.vector(res.dce.nolib$dce_pvalue),
    dce.tpm=as.vector(res.dce.tpm$dce_pvalue),
    dce.tpmlog=as.vector(res.dce.tpmlog$dce_pvalue),
    dce.lm=as.vector(res.dce.lm$dce_pvalue),
    rand=as.vector(res.rand$dce_pvalue),
    causaldag=as.vector(res.causaldag$dce_pvalue)
  )

  df.runtime <- data.frame(
    cor=time.cor,
    pcor=time.pcor,
    dce=time.dce,
    dce.log=time.dce.log,
    dce.latent=time.dce.latent,
    dce.lr=time.dce.lr,
    dce.nolib=time.dce.nolib,
    dce.tpm=time.dce.tpm,
    dce.tpmlog=time.dce.tpmlog,
    dce.lm=time.dce.lm,
    rand=time.rand,
    causaldag=time.causaldag
  )

  return(list(edges=df.res, pvalues=df.pvalues, runtime=df.runtime))
}
