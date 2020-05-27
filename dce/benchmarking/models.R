run.all.models <- function(
  wt.graph, wt.X,
  mt.graph, mt.X,
  wt.graph.perturbed,
  beta.magnitude,
  methods = c("cor", "pcor", "dce", "dce.lr", "rand"),
  lib.size = NULL
) {
  # for null model
  negweight.range <- c(-beta.magnitude, 0)
  posweight.range <- c(0, beta.magnitude)

  # ground truth
  ground.truth <- list(dce = trueEffects(mt.graph) - trueEffects(wt.graph))
  # reduce data for correlations
  graph <- as(wt.graph, "matrix")
  wt.X.cor <- wt.X[, seq_len(ncol(graph))]
  mt.X.cor <- mt.X[, seq_len(ncol(graph))]

  # correlations
  time.tmp <- Sys.time()
  if ("cor" %in% methods) {
    res.cor <- list(dce = cor(mt.X.cor) - cor(wt.X.cor))
    res.cor$dce.pvalue <- pcor_perm(wt.X.cor, mt.X.cor, fun = cor)
  } else {
    res.cor <- ground.truth
    res.cor$dce.pvalue <- ground.truth$dce*0
  }
  res.cor$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
  res.cor$dce.pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
  time.cor <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  time.tmp <- Sys.time()
  if ("pcor" %in% methods) {
    res.pcor <- list(dce = pcor(mt.X.cor) - pcor(wt.X.cor))
    res.pcor$dce.pvalue <- pcor_perm(wt.X.cor, mt.X.cor, fun = pcor)
  } else {
    res.pcor <- ground.truth
    res.pcor$dce.pvalue <- ground.truth$dce*0
  }
  res.pcor$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
  res.pcor$dce.pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
  time.pcor <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  # differential causal effects
  if (link.method == "log") {
    solver.args <- list(method = "glm.fit", link = "log")
  } else {
    solver.args <- list(method = "glm.dce.nb.fit", link = "identity")
  }

  time.tmp <- Sys.time()
  if ("dce" %in% methods) {
    res.dce <- dce::dce.nb(
      wt.graph.perturbed, wt.X, mt.X,
      adjustment.type = adjustment.type,
      solver.args = solver.args
    )
  } else {
    res.dce <- ground.truth
    res.dce$dce.pvalue <- ground.truth$dce*0
    res.dce$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
    res.dce$dce.pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
  }
  time.dce <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  time.tmp <- Sys.time()
  if ("dce.lr" %in% methods) {
    res.dce.lr <- dce::dce.nb(
      wt.graph.perturbed, wt.X, mt.X,
      adjustment.type = adjustment.type,
      solver.args = solver.args,#, test = "lr",
      lib.size = TRUE
    )
  } else {
    res.dce.lr <- ground.truth
    res.dce.lr$dce.pvalue <- ground.truth$dce*0
    res.dce.lr$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
    res.dce.lr$dce.pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
  }
  time.dce.lr <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  time.tmp <- Sys.time()
  if ("rand" %in% methods) {
    tmp <- as(wt.graph.perturbed, "matrix") * NA
    tmp[which(as(wt.graph.perturbed, "matrix") != 0)] = (
      runif(sum(as(wt.graph.perturbed, "matrix") != 0), negweight.range[1], posweight.range[2]) -
        runif(sum(as(wt.graph.perturbed, "matrix") != 0), negweight.range[1], posweight.range[2])
    )
    tmp.pvals <- as(wt.graph.perturbed, "matrix") * NA
    tmp.pvals[which(as(wt.graph.perturbed, "matrix") != 0)] <- runif(sum(as(wt.graph.perturbed, "matrix") != 0), 0, 1)
    res.rand <- list(dce = tmp, dce.pvalue = tmp.pvals)
  } else {
    res.rand <- ground.truth
    res.rand$dce.pvalue <- ground.truth$dce*0
    res.rand$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
    res.rand$dce.pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
  }
  time.rand <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  # gather results
  df.res <- data.frame(
    truth=as.vector(ground.truth$dce),
    cor=as.vector(res.cor$dce),
    pcor=as.vector(res.pcor$dce),
    dce=as.vector(res.dce$dce),
    dce.lr=as.vector(res.dce.lr$dce),
    rand=as.vector(res.rand$dce)
  )

  df.pvalues <- data.frame(
    truth=as.vector(ground.truth$dce),
    cor=as.vector(res.cor$dce.pvalue),
    pcor=as.vector(res.pcor$dce.pvalue),
    dce=as.vector(res.dce$dce.pvalue),
    dce.lr=as.vector(res.dce.lr$dce.pvalue),
    rand=as.vector(res.rand$dce.pvalue)
  )

  df.runtime <- data.frame(
    cor=time.cor,
    pcor=time.pcor,
    dce=time.dce,
    dce.lr=time.dce.lr,
    rand=time.rand
  )

  return(list(edges=df.res, pvalues=df.pvalues, runtime=df.runtime))
}
