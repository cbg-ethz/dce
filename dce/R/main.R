#' Compute the true casual effects of a simulated dag
#'
#' This function takes a DAG with edgeweights as input and computes
#' the causal effects of all nodes on all direct and indirect children in the
#' DAG. Alternatively see pcalg::causalEffect for pairwise computation.
#' @param g graphNEL object
#' @author Martin Pirkl
#' @return matrix of causal effects
#' @export
#' @importFrom pcalg causalEffect
#' @import graph tidyverse
#' @examples
trueEffects <- function(g) {
    n <- ncol(as(g, "matrix"))
    te <- matrix(0, n, n)
    for (i in seq_len(n-1)) {
        for (j in i:n) {
            te[i, j] <- causalEffect(g, j, i)
        }
    }
    return(te)
}

#' Compute the causal effects
#'
#' This function takes a DAG and gene expression data as input and computes
#' the causal effects of all nodes on all other nodes in the DAG.
#' @param graph DAG as a graphNEL object
#' @param df.expr genes expression profiles with rows as observations and
#' columns as variables (nodes)
#' @author Kim Jablonski
#' @return vector of causal effects
#' @export
#' @importFrom purrr map_dfc
#' @importFrom pcalg idaFast
#' @importFrom assertthat are_equal
#' @import graph tidyverse
#' @examples
#' dag <- matrix(c(0,0,0,1,0,0,0,1,0), 3)
#' colnames(dag) <- rownames(dag) <- seq_len(3)
#' dag <- as(dag, "graphNEL")
#' d <- matrix(rnorm(100*3), 100)
#' colnames(d) <- seq_len(3)
#' compute_causal_effects(dag, d)
compute_causal_effects <- function(graph, df.expr) {
  node.list <- nodes(graph)

  cov.mat <- cov(df.expr)
  assertthat::are_equal(dim(cov.mat)[[1]], length(node.list))

  n <- length(node.list)
  purrr::map_dfc(1:n, function(x) {
    pcalg::idaFast(x, 1:n, cov.mat, graph)
  })
}
#' Compute the differential causal effects
#'
#' This function takes a DAG and gene expression data of two distinct types
#' of populations (e.g. tumor versus normal) as input and computes
#' the differential causal effects of all nodes on all other nodes in the DAG
#' between the tow populations of samples.
#' @param graph.ctrl DAG as a graphNEL object of the control population
#' @param df.expr.ctrl genes expression profiles with rows as observations and
#' columns as variables (nodes) of the control population
#' @param graph.mut DAG as a graphNEL object of the tumor population
#' @param df.expr.mut genes expression profiles with rows as observations and
#' columns as variables (nodes) of the tumor population
#' @param either "full" for the full linear model using all sampels together to
#' directly compute the differential causal effect or "part" for two lienar
#' models on both conditions to first compute the causal effect and then the
#' difference
#' @param bootstrap if TRUE, uses a bootstrap or subsample routine to
#' compute the effects
#' @param runs bootstrap/subsampling runs
#' @param replace if TRUE, classic bootstrap, if FALSE subsampling without
#' replacement
#' @param frac the fraction of the data to sample, can either be one value
#' or a vector of two for different fractions for control and tumor,
#' respectively populations
#' @param bootMethod if "diff" (default) bottstraps the differential causal
#' effects, if "cause", bootstraps on the causal effects
#' @author Kim Jablonski & Martin Pirkl
#' @return vector of causal effects
#' @export
#' @importFrom purrr map_dfc
#' @importFrom pcalg idaFast
#' @importFrom assertthat are_equal
#' @import graph tidyverse
#' @examples
compute_differential_causal_effects <- function(graph.ctrl, df.expr.ctrl,
                                                graph.mut, df.expr.mut,
                                                method = "full",
                                                bootstrap = FALSE, runs = 100,
                                                replace = FALSE, frac = 0.5,
                                                bootMethod = 0, ...) {
    if (bootstrap) {
        if (method %in% "full") {
            strap <- 1
        }
        if (length(frac) == 1) { frac <- c(frac, frac) }
        if (ceiling(nrow(df.expr.ctrl)*frac[1]) < ncol(df.expr.ctrl)) {
            frac[1] <- ncol(df.expr.ctrl)/nrow(df.expr.ctrl)
            print("too few control samples for subsampling at this fraction; reset to:")
            print(frac[1])
        }
        if (ceiling(nrow(df.expr.mut)*frac[2]) < ncol(df.expr.mut)) {
            frac[2] <- ncol(df.expr.mut)/nrow(df.expr.mut)
            print("too few tumor samples for subsampling at this fraction; reset to:")
            print(frac[2])
        }
        if (ncol(df.expr.mut) == nrow(df.expr.mut)) {
            print("too few samples for sampling with replacement")
            replace <- FALSE
        }
        if (bootMethod %in% "diff") {
            dces <- 0
            for (b in seq_len(runs)) {
                if (!(method %in% "full")) {
                    df.expr.ctrl.sub <- df.expr.ctrl[
                        sample(seq_len(nrow(df.expr.ctrl)),
                               ceiling(nrow(df.expr.ctrl)*frac[1]),
                               replace = replace), ]
                    df.expr.mut.sub <- df.expr.mut[
                        sample(seq_len(nrow(df.expr.mut)),
                               ceiling(nrow(df.expr.mut)*frac[2]),
                               replace = replace), ]
                    dces <- dces +
                        compute_differential_causal_effects(graph.ctrl,
                                                            df.expr.ctrl.sub,
                                                            graph.mut,
                                                            df.expr.mut.sub)$dce
                } else {
                    dces <- dces + fulllin(graph.ctrl, df.expr.ctrl,
                                    graph.mut, df.expr.mut,
                                    ...)$dce
                }
            }
            res <- dces/runs
        } else {
            ce.ctrl <- ce.mut <- 0
            for (b in seq_len(runs)) {
                df.expr.ctrl.sub <- df.expr.ctrl[
                    sample(seq_len(nrow(df.expr.ctrl)),
                           ceiling(nrow(df.expr.ctrl)*frac[1]),
                           replace = replace), ]
                df.expr.mut.sub <- df.expr.mut[
                    sample(seq_len(nrow(df.expr.mut)),
                           ceiling(nrow(df.expr.mut)*frac[2]),
                           replace = replace), ]
                ce.ctrl <- ce.ctrl+compute_causal_effects(graph.ctrl,
                                                          df.expr.ctrl.sub)
                ce.mut <- ce.mut+compute_causal_effects(graph.mut,
                                                        df.expr.mut.sub)
            }
            ce.ctrl <- ce.ctrl/runs
            ce.mut <- ce.mut/runs
            res <- t(as.matrix(ce.ctrl - ce.mut))
        }
    } else {
        if (!(method %in% "full")) {
            ce.ctrl <- compute_causal_effects(graph.ctrl, df.expr.ctrl)
            ce.mut <- compute_causal_effects(graph.mut, df.expr.mut)
            res <- t(as.matrix(ce.ctrl - ce.mut))
        } else {
            res <- fulllin(graph.ctrl, df.expr.ctrl,
                                    graph.mut, df.expr.mut,
                                    ...)$dce
        }
    }
    if (!(method %in% "full")) {
        res <- list(dce = res, graph = graph.ctrl,
                    cen = t(as(ce.ctrl, "matrix")),
                    cet = t(as(ce.mut, "matrix")))
    } else {
        res <- list(dce = res, graph = graph.ctrl)
    }

    rownames(res$dce) <- nodes(graph.ctrl)
    colnames(res$dce) <- nodes(graph.ctrl)

    class(res) <- "dce"
    return(res)
}
#' Plot dce object
#'
#' This function takes a differential causal effects object and plots
#' the dag with the dces
#' @param x dce object
#' @param dec rounding to dec decimals
#' @author Martin Pirkl
#' @method plot dce
#' @return plot of dag and dces
#' @export
#' @examples
plot.dce <- function(x, dec=3, ...) {
    adj <- as(x$graph, "matrix")
    edge.weights <- x$dce[which(adj!=0)]

    e.max <- max(abs(edge.weights))
    color.mat <- colorRamp(c("red", "blue"))(scales::rescale(edge.weights, from=c(-e.max,e.max))) / 255

    # TODO: try using `mnem::plotDnf` to reduce "hairballing" of network
    GGally::ggnet2(
      adj,
      label=TRUE,
      edge.label=edge.weights %>% round(2) %>% as.character,
      edge.size=scales::rescale(abs(edge.weights), c(1, 3)),
      edge.color=apply(color.mat, 1, function(x) { rgb(x[1], x[2], x[3]) }),
      arrow.size=12, arrow.gap=0.025
    )
}
#' Simulation study
#'
#' This function takes several parameters to define a simulation study.
#' @param nodes numebr of genes
#' @param samples vector of length two for sample number of both conditions
#' @param runs simulations runs
#' @param mu mean expression
#' @param sd standard deviation
#' @param effRange vector of length four with lowest, second lowest, second highest
#' and highest possible effect the effects are uniformly drawn from
#' @param truePos the fraction of differential effects; truePos=1 means all
#' are potentially differential (see effRange); if truePos=0.5, 50% of
#' differential effects are set to 0.
#' @param perturb positive or negative frantion of edges to be added or
#' removed, respectively
#' @param corMeth method for the correlation accuracy (see ?cor)
#' @param prob edge probability; either probability or "runif" to draw a
#' probability in each run
#' @param verbose verbose output, if TRUE
#' @author Martin Pirkl
#' @return accuracy for several different methods and statistics
#' for the ground truth as a list of two arrays
#' @export
#' @importFrom nem transitive.reduction
#' @examples
simDce <- function(nodes=5, samples=c(10,10),runs=10,mu=0,sd=1,
                   effRange=c(-1,0,0,1),truePos=1,perturb=0,corMeth="p",
                   prob="runif",verbose=FALSE) {
    acc <- array(0, c(runs,5,2),
                 dimnames = list(runs = paste0("run_", seq_len(runs)),
                                 methods = c("dce", "random",
                                             "full linear",
                                             "simple correlation",
                                             "test"),
                                 metrics = c("correlation", "time")))
    gtnfeat <- array(0, c(runs, 4),
                     dimnames = list(runs = paste0("run_", seq_len(runs)),
                                     features = c("avg children/parents",
                                                  "max children",
                                                  "max parents",
                                                  "density")))
    n <- nodes
    m <- samples
    lB <- effRange[seq_len(2)]
    uB <- effRange[3:4]
    truepos <- truePos
    for (run in seq_len(runs)) {
        if (prob %in% "runif") {
            p2 <- runif(1)
        } else {
            p2 <- prob
        }
        normal <- randomDAG_2(n, p2, lB, uB)
        tumor <- newWeights(normal, lB, uB, truepos) # resample edge weights
        dn <- rmvDAG_2(m[2], normal, normpars = c(mu,sd))
        dt <- rmvDAG_2(m[1], tumor, normpars = c(mu,sd))
        gm <- as(normal, "matrix")
        gm[which(gm != 0)] <- 1
        cn <- trueEffects(normal)
        ct <- trueEffects(tumor)
        gtc <- nem::transitive.closure(gm, mat=TRUE)
        ## save features of gtn, which might correlate with accuracy:
        gtnfeat[run, 1] <- mean(apply(gm, 1, sum))
        gtnfeat[run, 2] <- max(apply(gm, 1, sum))
        gtnfeat[run, 3] <- max(apply(gm, 2, sum))
        gtnfeat[run, 4] <- sum(gm)
        ## ground truth:
        dcet <- (cn - ct)*gtc # gtn for differential causal effects
        dcegtn <- list(dce = dcet, graph = normal, dcefull = dcet)
        class(dcegtn) <- "dce"
        ## perturb network:
        if (perturb < 0) {
            adjn <- as(normal, "matrix")
            remedge <- sample(which(adjn != 0),
                              floor(sum(adjn != 0)*abs(perturb)))
            adjn[remedge] <- 0
            adjn[which(adjn != 0)] <- 1
            normal <- tumor <- as(adjn, "graphNEL")
        }
        if (perturb > 0) {
            adjn <- as(normal, "matrix")
            addedge <- sample(which(adjn == 0 & upper.tri(adjn)),
                              floor(sum(adjn == 0 & upper.tri(adjn))*perturb))
            adjn[addedge] <- 1
            adjn[which(adjn != 0)] <- 1
            normal <- tumor <- as(adjn, "graphNEL")
        }
        ## full linear model:
        start <- as.numeric(Sys.time())
        dcei <- compute_differential_causal_effects(
                   normal, dn,
                   tumor, dt, method = "full"
        )
        acc[run, 3, 2] <- as.numeric(Sys.time()) - start
        dceifl <- dcei
        dcei <- dcei$dce
        acc[run, 3, 1] <- cor(as.vector(dcet), as.vector(dcei), method = cormeth, use = "complete.obs")
        ## normal
        Cn <- cov(dn)
        Ct <- cov(dt)
        if (Matrix::rankMatrix(Cn) == nrow(Cn) & Matrix::rankMatrix(Ct) == nrow(Ct)) {
            start <- as.numeric(Sys.time())
            dcei <- compute_differential_causal_effects(
                normal, dn,
                tumor, dt, method = "normal"
            )
            acc[run, 1, 2] <- as.numeric(Sys.time()) - start
            dcein <- dcei
            dcei <- dcei$dce
            acc[run, 1, 1] <- cor(as.vector(dcet), as.vector(dcei), method = cormeth, use = "complete.obs")
        }
        ## simple correlation:
        start <- as.numeric(Sys.time())
        dcei <- (cor(dn) - cor(dt))*gtc
        acc[run, 4, 2] <- as.numeric(Sys.time()) - start
        dcec <- list(dce = dcei, graph = as(gtc, "graphNEL"), dcefull = dcei)
        dceic <- dcec
        class(dceic) <- "dce"
        dcec <- dcec$dce
        acc[run, 4, 1] <- cor(as.vector(dcet), as.vector(dcec), method = cormeth, use = "complete.obs")
        ## random base line:
        dcei <- dceicn <- dceict <- dcet
        start <- as.numeric(Sys.time())
        dcei[which(gtc != 0)] <- runif(sum(gtc != 0), lB[1], uB[2]) - runif(sum(gtc != 0), lB[1], uB[2])
        acc[run, 2, 2] <- as.numeric(Sys.time()) - start
        dcer <- list(dce = dcei, graph = as(gtc, "graphNEL"), dcefull = dcei)
        dceir <- dcer
        class(dceir) <- "dce"
        dcer <- dcer$dce
        coridx <- which(dcet != 0 | dcer != 0)
        acc[run, 2, 1] <- cor(as.vector(dcet), as.vector(dcer), method = cormeth, use = "complete.obs")
        if (verbose) {
            cat(paste0(run, "."))
        }
    }
    dceSim <- list(acc=acc,gtnFeat=gtnfeat)
    class(dceSim) <- "dceSim"
    return(dceSim)
}
#' Plot simulation study
#'
#' Takes a dceSim object and produces a figure.
#' @param x dceSim object
#' @param ... additional arguments for mnem::mnemBox
#' @author Martin Pirkl
#' @method plot dceSim
#' @return plot
#' @export
#' @importFrom nem transitive.reduction
#' @importFrom mnem mnemBox
#' @examples
plot.dceSim <- function(x, col = 1:4, showMeth = seq_len(4),
                        showFeat = 1, methNames = NULL, ...) {
    runs <- dim(x$acc)[1]
    if (is.null(methNames)) {
        methNames <- dimnames(x$acc)[[2]][showMeth]
    }
    par(mfrow=c(1,length(showFeat)))
    if (1 %in% showFeat) {
        myboxplot(x$acc[seq_len(runs), showMeth, 1], col = col,
                  main="Correlation", , xaxt = "n", ...)
        axis(1, seq_len(length(showMeth)), labels = methNames)
    }
    if (2 %in% showFeat) {
        myboxplot(x$acc[seq_len(runs), showMeth, 2], col = col,
                  main="Time", , xaxt = "n", ...)
        axis(1, seq_len(length(showMeth)), labels = methNames)
    }
    if (3 %in% showFeat) {
        myboxplot(x$gtnFeat[seq_len(runs), seq_len(4)], col = col,
                  main="Ground truth Features", , xaxt = "n", ...)
        axis(1, seq_len(length(showMeth)),
             labels = dimnames(x$gtnFeat)[[2]][seq_len(4)])
    }
}
