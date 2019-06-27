library(pcalg)
library(mnem)
library(graph)

## add noise options to the pcalg data generating function:

rmvDAG_2 <- function (n, dag, errDist = c("normal", "cauchy", "t4", "mix",
    "mixt3", "mixN100"), normpars = c(0,1), mix = 0.1, errMat = NULL, back.compatible = FALSE,
    use.node.names = !back.compatible) {
    stopifnot(is(dag, "graph"), (p <- length(nodes(dag))) >=
        2)
    weightMatrix <- if (back.compatible)
        wgtMatrix.0(dag)
    else wgtMatrix(dag)
    nonZeros <- which(weightMatrix != 0, arr.ind = TRUE)
    if (nrow(nonZeros) > 0) {
        if (any(nonZeros[, 1] - nonZeros[, 2] < 0) || any(diag(weightMatrix) !=
            0))
            stop("Input DAG must be topologically ordered!")
    }
    errDist <- match.arg(errDist)
    if (grepl("^mix", errDist))
        eMat <- function(outs) {
            X <- c(rnorm(n * p - length(outs)), outs)
            matrix(sample(X), nrow = n)
        }
    if (is.null(errMat)) {
        errMat <- switch(errDist, normal = matrix(rnorm(n * p, normpars[1], normpars[2]),
            nrow = n), cauchy = matrix(rcauchy(n * p), nrow = n),
            t4 = matrix(rt(n * p, df = 4), nrow = n), mix = eMat(rcauchy(round(mix *
                n * p))), mixt3 = eMat(rt(round(mix * n * p),
                df = 3)), mixN100 = eMat(rnorm(round(mix * n *
                p), sd = 10)))
    }
    else {
        stopifnot(!is.null(dim.eM <- dim(errMat)), dim.eM ==
            c(n, p), is.numeric(errMat))
    }
    if (use.node.names)
        colnames(errMat) <- nodes(dag)
    if (sum(weightMatrix) > 0) {
        X <- errMat
        for (j in 2:p) {
            ij <- 1:(j - 1)
            X[, j] <- X[, j] + X[, ij, drop = FALSE] %*% weightMatrix[j,
                ij]
        }
        X
    }
    else errMat
}

## put noise (here resample) on the weights of a dag (for tumor samples):

newWeights <- function(g, lB = -1, uB = 1, tp = 0.5) {
    n <- length(g@nodes)
    w <- g@edgeData@data
    arediff <- sample(seq_len(length(w)), floor(tp*length(w)))
    for (i in arediff) {
        w[[i]]$weight <- runif(1, lB, uB)
    }
    g@edgeData@data <- w
    return(g)
}
