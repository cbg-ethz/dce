compute.tpm <- function(counts, len.kb = 1) {
  # TODO: is len.kb == 1 reasonable in our case?
  rpk <- counts / len.kb
  pm.scale <- rowSums(rpk) / 1e6
  return(rpk / pm.scale)
}


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
  wt.X.cor <- compute.tpm(wt.X)[, seq_len(ncol(graph))]
  mt.X.cor <- compute.tpm(mt.X)[, seq_len(ncol(graph))]

  # correlations
  time.tmp <- Sys.time()
  if (is.null(methods) || "cor" %in% methods) {
    res.cor <- list(dce = cor(mt.X.cor) - cor(wt.X.cor))
    res.cor$dce_pvalue <- permutation_test(wt.X.cor, mt.X.cor, fun = cor)
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
    res.pcor$dce_pvalue <- permutation_test(wt.X.cor, mt.X.cor, fun = pcor)
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
    res.dce.tpm <- dce::dce_nb(
      wt.graph.perturbed, compute.tpm(wt.X), compute.tpm(mt.X),
      adjustment_type = adjustment.type,
      effect_type = effect.type,
      solver_args = solver.args,
      lib_size = FALSE,
      deconfounding = latent
    )
  } else {
    res.dce.tpm <- ground.truth
    res.dce.tpm$dce_pvalue <- ground.truth$dce*0
    res.dce.tpm$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
    res.dce.tpm$dce_pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
  }
  time.dce.tpm <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))
  
  time.tmp <- Sys.time()
  if (is.null(methods) || "dce.tpm.nolatent" %in% methods) {
    res.dce.tpm.nolatent <- dce::dce_nb(
      wt.graph.perturbed, compute.tpm(wt.X), compute.tpm(mt.X),
      adjustment_type = adjustment.type,
      effect_type = effect.type,
      solver_args = solver.args,
      lib_size = FALSE
    )
  } else {
    res.dce.tpm.nolatent <- ground.truth
    res.dce.tpm.nolatent$dce_pvalue <- ground.truth$dce*0
    res.dce.tpm.nolatent$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
    res.dce.tpm.nolatent$dce_pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
  }
  time.dce.tpm.nolatent <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))
  
  time.tmp <- Sys.time()
  if (is.null(methods) || "dce.lm.tpm" %in% methods) {
    res.dce.lm.tpm <- dce(
      wt.graph.perturbed, compute.tpm(wt.X), compute.tpm(mt.X),
      solver = "lm",
      adjustment_type = adjustment.type,
      effect_type = effect.type,
      lib_size = FALSE,
      deconfounding = latent
    )
  } else {
    res.dce.lm.tpm <- ground.truth
    res.dce.lm.tpm$dce_pvalue <- ground.truth$dce*0
    res.dce.lm.tpm$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
    res.dce.lm.tpm$dce_pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
  }
  time.dce.lm.tpm <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  # LDGM
  time.tmp <- Sys.time()
  if (is.null(methods) || "ldgm" %in% methods) {
    res.ldgm <- list(dce = LDGM(log(wt.X.cor+1),log(mt.X.cor+1)))
    res.ldgm$dce_pvalue <- permutation_test(log(wt.X.cor+1),log(mt.X.cor+1),fun=LDGM)
  } else {
    res.ldgm <- ground.truth
    res.ldgm$dce_pvalue <- ground.truth$dce*0
    res.ldgm$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
    res.ldgm$dce_pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
  }
  time.ldgm <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  # Diff with FastGGM
  time.tmp <- Sys.time()
  if (is.null(methods) || "fggm" %in% methods) {
    res.fggm <- FastGGM_Diff(log(wt.X.cor+1),log(mt.X.cor+1))
  } else {
    res.fggm <- ground.truth
    res.fggm$dce_pvalue <- ground.truth$dce*0
    res.fggm$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
    res.fggm$dce_pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
  }
  time.fggm <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

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
    dname.tmp <- paste0("tmp.",runif(1),".causaldag/")

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
    dce.nolib=as.vector(res.dce.nolib$dce),
    dce.tpm=as.vector(res.dce.tpm$dce),
    dce.tpm.nolatent=as.vector(res.dce.tpm.nolatent$dce),
    dce.lm.tpm=as.vector(res.dce.lm.tpm$dce),
    ldgm=as.vector(res.ldgm$dce),
    fggm=as.vector(res.fggm$dce),
    rand=as.vector(res.rand$dce),
    causaldag=as.vector(res.causaldag$dce)
  )

  df.pvalues <- data.frame(
    truth=as.vector(ground.truth$dce),
    cor=as.vector(res.cor$dce_pvalue),
    pcor=as.vector(res.pcor$dce_pvalue),
    dce=as.vector(res.dce$dce_pvalue),
    dce.nolib=as.vector(res.dce.nolib$dce_pvalue),
    dce.tpm=as.vector(res.dce.tpm$dce_pvalue),
    dce.tpm.nolatent=as.vector(res.dce.tpm.nolatent$dce_pvalue),
    dce.lm.tpm=as.vector(res.dce.lm.tpm$dce_pvalue),
    ldgm=as.vector(res.ldgm$dce_pvalue),
    fggm=as.vector(res.fggm$dce_pvalue),
    rand=as.vector(res.rand$dce_pvalue),
    causaldag=as.vector(res.causaldag$dce_pvalue)
  )

  df.runtime <- data.frame(
    cor=time.cor,
    pcor=time.pcor,
    dce=time.dce,
    dce.nolib=time.dce.nolib,
    dce.tpm=time.dce.tpm,
    dce.tpm.nolatent=time.dce.tpm.nolatent,
    dce.lm.tpm=time.dce.lm.tpm,
    ldgm=time.ldgm,
    fggm=time.fggm,
    rand=time.rand,
    causaldag=time.causaldag
  )

  return(list(edges=df.res, pvalues=df.pvalues, runtime=df.runtime))
}
