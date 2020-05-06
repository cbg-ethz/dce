run.all.models <- function(
  wt.graph, wt.X,
  mt.graph, mt.X,
  wt.graph.perturbed,
  beta.magnitude
) {
  # for null model
  negweight.range <- c(-beta.magnitude, 0)
  posweight.range <- c(0, beta.magnitude)

  # ground truth
  ground.truth <- list(dce = trueEffects(mt.graph) - trueEffects(wt.graph))

  # correlations
  time.tmp <- Sys.time()
  res.cor <- list(dce = cor(mt.X) - cor(wt.X))
  res.cor$dce.pvalue <- pcor_perm(wt.X, mt.X, fun = cor)
  time.cor <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  time.tmp <- Sys.time()
  res.pcor <- list(dce = pcor(mt.X) - pcor(wt.X))
  res.pcor$dce.pvalue <- pcor_perm(wt.X, mt.X, fun = pcor)
  time.pcor <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  # differential causal effects
  if (link.method == "log") {
    solver.args <- list(method = "glm.fit", link = "log")
  } else {
    solver.args <- list(method = "glm.dce.nb.fit", link = "identity")
  }

  time.tmp <- Sys.time()
  res.dce <- dce::dce.nb(
    wt.graph.perturbed, wt.X, mt.X,
    adjustment.type = adjustment.type,
    solver.args = solver.args
  )
  time.dce <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  time.tmp <- Sys.time()
  res.dce.lr <- dce::dce.nb(
    wt.graph.perturbed, wt.X, mt.X,
    adjustment.type = adjustment.type,
    solver.args = solver.args, test = "lr"
  )
  time.dce.lr <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  time.tmp <- Sys.time()
  tmp <- as.matrix(ground.truth$dce)
  tmp[which(as.matrix(ground.truth$dce) != 0)] = (
    runif(sum(as.matrix(ground.truth$dce) != 0), negweight.range[1], posweight.range[2]) -
      runif(sum(as.matrix(ground.truth$dce) != 0), negweight.range[1], posweight.range[2])
  )
  tmp.pvals <- tmp*0
  tmp.pvals[which(as.matrix(ground.truth$dce) != 0)] <- runif(sum(as.matrix(ground.truth$dce) != 0), 0, 1)
  res.rand <- list(dce = tmp, dce.pvalue = tmp.pvals)
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
