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
  graph1 <- as_adjmat(wt.graph)
  wt.X.cor <- compute.tpm(wt.X)[, seq_len(ncol(graph))]
  mt.X.cor <- compute.tpm(mt.X)[, seq_len(ncol(graph))]

  # correlations
  time.tmp <- Sys.time()
  if (is.null(methods) || "cor" %in% methods) {
    res.cor <- list(dce = cor(mt.X.cor, method = 'spearman') - cor(wt.X.cor, method = 'spearman'))
    res.cor$dce_pvalue <- permutation_test(wt.X.cor, mt.X.cor, fun = cor, method = 'spearman')
  } else {
    res.cor <- ground.truth
    res.cor$dce_pvalue <- ground.truth$dce*0
  }
  time.cor <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  time.tmp <- Sys.time()
  if (is.null(methods) || "pcor" %in% methods) {
    res.pcor <- list(dce = pcor(mt.X.cor, method = 'spearman') - pcor(wt.X.cor, method = 'spearman'))
    res.pcor$dce_pvalue <- permutation_test(wt.X.cor, mt.X.cor, fun = pcor, method = 'spearman')
  } else {
    res.pcor <- ground.truth
    res.pcor$dce_pvalue <- ground.truth$dce*0
  }
  time.pcor <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  time.tmp <- Sys.time()
  if (is.null(methods) || "pcorz" %in% methods) {
    res.pcorz <- list(dce = pcor(mt.X.cor, g = graph1, method = 'spearman') - pcor(wt.X.cor, g = graph1, method = 'spearman', adjustment_type = adjustment.type))
    res.pcorz$dce_pvalue <- permutation_test(wt.X.cor, mt.X.cor, fun = pcor, method = 'spearman', g = graph1, adjustment_type = adjustment.type)
  } else {
    res.pcorz <- ground.truth
    res.pcorz$dce_pvalue <- ground.truth$dce*0
  }
  time.pcorz <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

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
  }
  time.dce <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  time.tmp <- Sys.time()
  if (is.null(methods) || "dce.nolib" %in% methods) {
    res.dce.nolib <- dce(
      wt.graph.perturbed, wt.X, mt.X,
      solver = "lm",
      adjustment_type = adjustment.type,
      effect_type = effect.type,
      lib_size = FALSE
    )
  } else {
    res.dce.nolib <- ground.truth
    res.dce.nolib$dce_pvalue <- ground.truth$dce*0
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
  }
  time.dce.lm.tpm <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  time.tmp <- Sys.time()
  if (is.null(methods) || "dce.lm.tpm.HC" %in% methods) {
    res.dce.lm.tpm.HC <- dce(
      wt.graph.perturbed, compute.tpm(wt.X), compute.tpm(mt.X),
      solver = "lm",
      adjustment_type = adjustment.type,
      effect_type = effect.type,
      lib_size = FALSE,
      deconfounding = latent,
      test = 'vcovHC'
    )
  } else {
    res.dce.lm.tpm.HC <- ground.truth
    res.dce.lm.tpm.HC$dce_pvalue <- ground.truth$dce*0
  }
  time.dce.lm.tpm.HC <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  time.tmp <- Sys.time()
  if (is.null(methods) || "dce.lm.tpm.nolatent" %in% methods) {
    res.dce.lm.tpm.nolatent <- dce(
      wt.graph.perturbed, compute.tpm(wt.X), compute.tpm(mt.X),
      solver = "lm",
      adjustment_type = adjustment.type,
      effect_type = effect.type,
      lib_size = FALSE
    )
  } else {
    res.dce.lm.tpm.nolatent <- ground.truth
    res.dce.lm.tpm.nolatent$dce_pvalue <- ground.truth$dce*0
  }
  time.dce.lm.tpm.nolatent <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  # LDGM
  time.tmp <- Sys.time()
  if (is.null(methods) || "ldgm" %in% methods) {
    res.ldgm <- list(dce = LDGM(log(wt.X.cor+1),log(mt.X.cor+1)))
    res.ldgm$dce_pvalue <- permutation_test(log(wt.X.cor+1),log(mt.X.cor+1),fun=LDGM)
  } else {
    res.ldgm <- ground.truth
    res.ldgm$dce_pvalue <- ground.truth$dce*0
  }
  time.ldgm <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  # Diff with FastGGM
  time.tmp <- Sys.time()
  if (is.null(methods) || "fggm" %in% methods) {
    res.fggm <- FastGGM_Diff(log(wt.X.cor+1),log(mt.X.cor+1))
  } else {
    res.fggm <- ground.truth
    res.fggm$dce_pvalue <- ground.truth$dce*0
  }
  time.fggm <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  # null models
  time.tmp <- Sys.time()
  if (is.null(methods) || "rand" %in% methods) {
    tmp.graph.perturbed <- as(wt.graph.perturbed, "matrix")
    tmp.graph.perturbed <- tmp.graph.perturbed[naturalorder(rownames(tmp.graph.perturbed)), naturalorder(colnames(tmp.graph.perturbed))]
    tmp <- tmp.graph.perturbed * NA
    tmp[which(tmp.graph.perturbed != 0)] = (
      runif(sum(tmp.graph.perturbed != 0), negweight.range[1], posweight.range[2]) -
        runif(sum(tmp.graph.perturbed != 0), negweight.range[1], posweight.range[2])
    )
    tmp.pvals <- tmp.graph.perturbed * NA
    tmp.pvals[which(tmp.graph.perturbed != 0)] <- runif(sum(tmp.graph.perturbed != 0), 0, 1)
    res.rand <- list(dce = tmp, dce_pvalue = tmp.pvals)
  } else {
    res.rand <- ground.truth
    res.rand$dce_pvalue <- ground.truth$dce*0
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
    }
  } else {
    res.causaldag <- ground.truth
    res.causaldag$dce_pvalue <- ground.truth$dce*0
  }
  time.causaldag <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  # gather results
  df.res <- data.frame(
    truth=as.vector(ground.truth$dce),
    cor=as.vector(res.cor$dce),
    pcor=as.vector(res.pcor$dce),
    pcorz=as.vector(res.pcorz$dce),
    dce=as.vector(res.dce$dce),
    dce.nolib=as.vector(res.dce.nolib$dce),
    dce.tpm=as.vector(res.dce.tpm$dce),
    dce.tpm.nolatent=as.vector(res.dce.tpm.nolatent$dce),
    dce.lm.tpm=as.vector(res.dce.lm.tpm$dce),
    dce.lm.tpm.HC=as.vector(res.dce.lm.tpm.HC$dce),
    dce.lm.tpm.nolatent=as.vector(res.dce.lm.tpm.nolatent$dce),
    ldgm=as.vector(res.ldgm$dce),
    fggm=as.vector(res.fggm$dce),
    rand=as.vector(res.rand$dce),
    causaldag=as.vector(res.causaldag$dce)
  )

  df.pvalues <- data.frame(
    truth=as.vector(ground.truth$dce),
    cor=as.vector(res.cor$dce_pvalue),
    pcor=as.vector(res.pcor$dce_pvalue),
    pcorz=as.vector(res.pcorz$dce_pvalue),
    dce=as.vector(res.dce$dce_pvalue),
    dce.nolib=as.vector(res.dce.nolib$dce_pvalue),
    dce.tpm=as.vector(res.dce.tpm$dce_pvalue),
    dce.tpm.nolatent=as.vector(res.dce.tpm.nolatent$dce_pvalue),
    dce.lm.tpm=as.vector(res.dce.lm.tpm$dce_pvalue),
    dce.lm.tpm.HC=as.vector(res.dce.lm.tpm.HC$dce_pvalue),
    dce.lm.tpm.nolatent=as.vector(res.dce.lm.tpm.nolatent$dce_pvalue),
    ldgm=as.vector(res.ldgm$dce_pvalue),
    fggm=as.vector(res.fggm$dce_pvalue),
    rand=as.vector(res.rand$dce_pvalue),
    causaldag=as.vector(res.causaldag$dce_pvalue)
  )

  df.runtime <- data.frame(
    cor=time.cor,
    pcor=time.pcor,
    pcorz=time.pcorz,
    dce=time.dce,
    dce.nolib=time.dce.nolib,
    dce.tpm=time.dce.tpm,
    dce.tpm.nolatent=time.dce.tpm.nolatent,
    dce.lm.tpm=time.dce.lm.tpm,
    dce.lm.tpm.HC=time.dce.lm.tpm.HC,
    dce.lm.tpm.nolatent=time.dce.lm.tpm.nolatent,
    ldgm=time.ldgm,
    fggm=time.fggm,
    rand=time.rand,
    causaldag=time.causaldag
  )

  return(list(edges=df.res, pvalues=df.pvalues, runtime=df.runtime))
}
