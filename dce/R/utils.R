library(pcalg)
library(mnem)
library(graph)

Gsolve <- function(a, b, ...) {
    x <- t(solve(diag(1, nrow(a)), Ginv(a)%*%b))
    return(x)
}

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

## stronger effects
randomDAG_2 <- function (n, prob, lB = 0.1, uB = 1, V = as.character(1:n))
{
    stopifnot(n >= 2, is.numeric(prob), length(prob) == 1, 0 <=
        prob, prob <= 1, is.numeric(lB), is.numeric(uB), lB <=
                                                         uB)
    if (length(uB) == 1) { uB <- c(0, uB) }
    if (length(lB) == 1) { lB <- c(lB, 0) }
    edL <- vector("list", n)
    nmbEdges <- 0L
    for (i in seq_len(n - 2)) {
        listSize <- rbinom(1, n - i, prob)
        nmbEdges <- nmbEdges + listSize
        edgeList <- sample(seq(i + 1, n), size = listSize)
        posedges <- sample(seq_len(length(edgeList)), ceiling(length(edgeList)/2))
        negedges <- seq_len(length(edgeList))[-posedges]
        weightList <- numeric(length(edgeList))
        weightList[posedges] <- runif(length(posedges), min = uB[1], max = uB[2])
        if (length(negedges) > 0) {
            weightList[negedges] <- runif(length(negedges), min = lB[1], max = lB[2])
        }
        edL[[i]] <- list(edges = edgeList, weights = weightList)
    }
    listSize <- rbinom(1, 1, prob)
    if (listSize > 0) {
        nmbEdges <- nmbEdges + 1
        edgeList <- n
        weightList <- runif(1, min = lB, max = uB)
    }
    else {
        edgeList <- integer(0)
        weightList <- numeric(0)
    }
    edL[[n - 1]] <- list(edges = edgeList, weights = weightList)
    if (nmbEdges > 0) {
        edL[[n]] <- list(edges = integer(0), weights = numeric(0))
        names(edL) <- V
        new("graphNEL", nodes = V, edgeL = edL, edgemode = "directed")
    }
    else new("graphNEL", nodes = V, edgemode = "directed")
}

## put noise (here resample) on the weights of a dag (for tumor samples):

newWeights <- function(g, lB = -1, uB = 1, tp = 0.5) {
    if (length(uB) == 1) { uB <- c(0, uB) }
    if (length(lB) == 1) { lB <- c(lB, 0) }
    n <- length(g@nodes)
    w <- g@edgeData@data
    arediff <- sample(seq_len(length(w)), floor(tp*length(w)))
    for (i in arediff) {
        die <- sample(c(0,1), 1)
        if (die) {
            w[[i]]$weight <- runif(1, lB[1], lB[2])
        } else {
            w[[i]]$weight <- runif(1, uB[1], uB[2])
        }
    }
    g@edgeData@data <- w
    return(g)
}
