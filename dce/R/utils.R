#' @importFrom matlib Ginv
#' @noRd
Gsolve <- function(a, b=NULL, ...) {
    if (!is.null(b)) {
        x <- t(solve(diag(1, nrow(a)), matlib::Ginv(a, ...)%*%b))
    } else {
        x <- matlib::Ginv(a, ...)
    }
    return(x)
}
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
fulllin <- function(g1, d1, g2, d2, conf = TRUE, diff = 1, ...) {
    mat1 <- as(g1, "matrix")
    mat2 <- as(g2, "matrix")
    mat1[which(mat1 != 0)] <- 1
    mat2[which(mat2 != 0)] <- 1
    dagtc <- mnem:::mytc(mat1)
    df <- rbind(d1, d2)
    colnames(df) <- paste0("X", seq_len(ncol(df)))
    df <- as.data.frame(cbind(df,
                              N = c(rep(1, nrow(d1)),
                                    rep(0, nrow(d2)))))

    n <- length(nodes(g1))

    if (diff) {
        dce <- mat1*0
        for (i in seq_len(n)) {
            for (j in seq_len(n)) {
                if (dagtc[i, j] == 1 & i != j) {
                    Z <- pcalg::backdoor(mat1, i, j, type = "dag")
                    NX <- df[, i] * df$N
                    NZ <- df[, Z] * df$N
                    X <- df[, i]
                    Y <- df[, j]
                    Z <- df[, Z]
                    if (length(Z) > 0 & conf) {
                        C <- cov(cbind(Y, NX, NZ, X, Z))
                    } else {
                        C <- cov(cbind(Y, NX, X))
                    }
                    if (Matrix::rankMatrix(C) < nrow(C)) {
                        betas <- Gsolve(C[2:nrow(C), 2:ncol(C)], C[2:nrow(C), 1])
                    } else {
                        betas <- solve(C[2:nrow(C), 2:ncol(C)], C[2:nrow(C), 1])
                    }
                    dce[i, j] <- betas[1]
                }
            }
        }
    }
    gtc <- as(g1, "matrix")
    gtc[which(gtc != 0)] <- 1
    gtc <- mnem:::mytc(gtc)
    diag(gtc) <- 0
    res <- list(dce = dce*gtc, graph = g1, dcefull = dce)
    class(res) <- "dce"
    return(res)
}
