run.all.models <- function(
  wt.graph, wt.X,
  mt.graph, mt.X,
  wt.graph.perturbed,
  beta.magnitude,
  methods = c("cor", "pcor", "dce", "dce.lr", "rand", "causaldag"),
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
      solver.args = solver.args,
      lib.size = TRUE#, conservative = TRUE
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
      solver.args = solver.args, test = "lr",
      lib.size = TRUE,
      latent = 1
    )
  } else {
    res.dce.lr <- ground.truth
    res.dce.lr$dce.pvalue <- ground.truth$dce*0
    res.dce.lr$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
    res.dce.lr$dce.pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
  }
  time.dce.lr <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  # null models
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

  # other methods
  time.tmp <- Sys.time()
  if ("causaldag" %in% methods) {
    # run causaldag
    dname.tmp <- "tmp.causaldag/"

    unlink(dname.tmp, recursive = TRUE)
    dir.create(dname.tmp, recursive = TRUE)

    write.csv(wt.X[, nodes(wt.graph.perturbed)], file.path(dname.tmp, "X_wt.csv"), quote=FALSE)
    write.csv(wt.X[, nodes(wt.graph.perturbed)], file.path(dname.tmp, "X_mt.csv"), quote=FALSE)

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
        dce.pvalue = as.matrix(df.causaldag)
      )

      res.causaldag$dce.pvalue[res.causaldag$dce.pvalue != 0] <- 0.01 # something significant...
      res.causaldag$dce.pvalue[res.causaldag$dce.pvalue == 0] <- 0.99

      res.causaldag$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
      res.causaldag$dce.pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
    }
  } else {
    res.causaldag <- ground.truth
    res.causaldag$dce.pvalue <- ground.truth$dce*0
    res.causaldag$dce[as(wt.graph.perturbed, "matrix") == 0] <- NA
    res.causaldag$dce.pvalue[as(wt.graph.perturbed, "matrix") == 0] <- NA
  }
  time.causaldag <- as.integer(difftime(Sys.time(), time.tmp, units = "secs"))

  # gather results
  df.res <- data.frame(
    truth=as.vector(ground.truth$dce),
    cor=as.vector(res.cor$dce),
    pcor=as.vector(res.pcor$dce),
    dce=as.vector(res.dce$dce),
    dce.lr=as.vector(res.dce.lr$dce),
    rand=as.vector(res.rand$dce),
    causaldag=as.vector(res.causaldag$dce)
  )

  df.pvalues <- data.frame(
    truth=as.vector(ground.truth$dce),
    cor=as.vector(res.cor$dce.pvalue),
    pcor=as.vector(res.pcor$dce.pvalue),
    dce=as.vector(res.dce$dce.pvalue),
    dce.lr=as.vector(res.dce.lr$dce.pvalue),
    rand=as.vector(res.rand$dce.pvalue),
    causaldag=as.vector(res.causaldag$dce.pvalue)
  )

  df.runtime <- data.frame(
    cor=time.cor,
    pcor=time.pcor,
    dce=time.dce,
    dce.lr=time.dce.lr,
    rand=time.rand,
    causaldag=time.causaldag
  )

  return(list(edges=df.res, pvalues=df.pvalues, runtime=df.runtime))
}
