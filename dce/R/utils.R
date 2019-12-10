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
#' Compute true positive/... counts
#'
#' Useful for performance evaluations
#' @param truth Ground truth
#' @param inferred Computed results
#' @param cutoff Threshold for classification
#' @author Hans Wurst
#' @return data.frame
#' @export
#' @examples
#' get_prediction_counts(c(1,0), c(1,1))
get_prediction_counts <- function(truth, inferred, cutoff = 0.5) {
    tp <- sum(
        abs(truth) > cutoff &
        abs(inferred) > cutoff &
        sign(truth) == sign(inferred) & inferred != 0
    )
    fn <- sum(
        abs(truth) > cutoff &
        abs(inferred) <= cutoff &
        inferred != 0
    )
    fp <- sum(
        abs(truth) <= cutoff &
        abs(inferred) > cutoff |
        (
            abs(truth) > cutoff &
            abs(inferred) > cutoff &
            sign(truth) != sign(inferred)
        ) &
        inferred != 0
    )
    tn <- sum(
        abs(truth) <= cutoff &
        abs(inferred) <= cutoff &
        inferred != 0
    )

    return(data.frame(
        true.positive=tp,
        false.positive=fp,
        true.negative=tn,
        false.negative=fn
    ))
}
#' Simulate data
#'
#' Generate data for given DAG.
#' @param dag Graph to simulate on
#' @param n Number of samples
#' @param errDist Noise distribution
#' @param normpars Noise parameters
#' @param mix Another noise parameter
#' @param errMat Initial error matrix
#' @param back.compatible Legacy compatibility
#' @param use.node.names Keep node names
#' @author Martin Pirkl
#' @return graph
#' @export
#' @importFrom methods is
#' @importFrom pcalg wgtMatrix
#' @importFrom MASS rnegbin
#' @examples
#' dag <- create_random_DAG(30, 0.2)
#' X <- simulate_data(dag)
simulate_data <- function (
    dag, n = 100,
    errDist = c("normal", "cauchy", "t4", "mix",
    "mixt3", "mixN100", "nbinom"), normpars = c(0,1),
    mix = 0.1,
    errMat = NULL, back.compatible = FALSE,
    use.node.names = !back.compatible
) {
    stopifnot(is(dag, "graph"), (p <- length(nodes(dag))) >= 2)
    # TODO: if(back.compatible) wgtMatrix.0(dag)
    weightMatrix <- pcalg::wgtMatrix(dag)
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
        errMat <- switch(
            errDist,
            normal = matrix(rnorm(n * p, normpars[1], normpars[2]), nrow = n),
            nbinom = matrix(rnegbin(n * p, mu = normpars[1],
                                    theta = normpars[2]), nrow = n),
            cauchy = matrix(rcauchy(n * p), nrow = n),
            t4 = matrix(rt(n * p, df = 4), nrow = n),
            mix = eMat(rcauchy(round(mix * n * p))),
            mixt3 = eMat(rt(round(mix * n * p), df = 3)),
            mixN100 = eMat(rnorm(round(mix * n * p), sd = 10))
        )
    }
    else {
        stopifnot(!is.null(dim.eM <- dim(errMat)), dim.eM ==
            c(n, p), is.numeric(errMat))
    }
    if (use.node.names) colnames(errMat) <- nodes(dag)
    if (sum(weightMatrix) > 0) {
        X <- errMat
        for (j in 2:p) {
            ij <- seq_len(j - 1)
            if (errDist %in% "nbinom") {
                mu <- X[, j] + X[, ij, drop = FALSE] %*% weightMatrix[j,
                    ij]
                X[, j] <- rnegbin(n, mu = mu, theta = normpars[2])
                X[is.na(X[, j]), j] <- 0
            } else {
                X[, j] <- X[, j] + X[, ij, drop = FALSE] %*% weightMatrix[j,
                    ij]
            }
        }
        ## add drop ins (should only affect nbinom)
        Genes0 <- which(apply(X, 2, sum) == 0)
        X[, Genes0] <- rnegbin(length(Genes0), mu = 1, theta = normpars[2])
        return(X)
    }
    else errMat
}
#' Create random DAG (topologically ordered)
#'
#' Creates a DAG according to given parameters.
#' @param n Number of nodes
#' @param prob Probability of creating an edge
#' @param lB Lower bound for edge weights
#' @param uB Upper bound for edge weights
#' @param node.labels Node labels
#' @author Martin Pirkl
#' @return graph
#' @export
#' @importFrom methods new
#' @examples
#' dag <- create_random_DAG(30, 0.2)
create_random_DAG <- function (
    n, prob,
    lB = -1, uB = 1,
    node.labels = as.character(seq_len(n))
) {
    stopifnot(
        n >= 2, is.numeric(prob), length(prob) == 1, 0 <= prob, prob <= 1,
        is.numeric(lB), is.numeric(uB), lB <= uB
    )
    if (length(uB) == 1) { uB <- c(0, uB) }
    if (length(lB) == 1) { lB <- c(lB, 0) }
    if (lB[1] - lB[2] == 0) {
        lB <- uB
    }
    if (uB[1] - uB[2] == 0) {
        uB <- lB
    }
    edL <- vector("list", n)
    nmbEdges <- 0L
    for (i in seq_len(n - 2)) {
        listSize <- rbinom(1, n - i, prob)
        nmbEdges <- nmbEdges + listSize
        edgeList <- sample(seq(i + 1, n), size = listSize)
        posedges <- sample(
            seq_len(length(edgeList)), ceiling(length(edgeList)/2)
        )
        negedges <- seq_len(length(edgeList))[-posedges]
        weightList <- numeric(length(edgeList))
        weightList[posedges] <- runif(
            length(posedges), min = uB[1], max = uB[2]
        )
        if (length(negedges) > 0) {
            weightList[negedges] <- runif(
                length(negedges), min = lB[1], max = lB[2]
            )
        }
        edL[[i]] <- list(edges = edgeList, weights = weightList)
    }
    listSize <- rbinom(1, 1, prob)
    if (listSize > 0) {
        nmbEdges <- nmbEdges + 1
        edgeList <- n
        weightList <- runif(1, min = lB[1], max = uB[2])
    }
    else {
        edgeList <- integer(0)
        weightList <- numeric(0)
    }
    edL[[n - 1]] <- list(edges = edgeList, weights = weightList)
    if (nmbEdges > 0) {
        edL[[n]] <- list(edges = integer(0), weights = numeric(0))
        names(edL) <- node.labels
        new("graphNEL", nodes = node.labels, edgeL = edL, edgemode = "directed")
    }
    else new("graphNEL", nodes = node.labels, edgemode = "directed")
}
#' Resample network edge weights
#'
#' Takes a graph and modifies edge weights.
#' @param g original graph
#' @param lB lower bound for new edge weights
#' @param uB upper bound for new edge weights
#' @param tp fraction of edge weights which will be modified
#' @author Martin Pirkl
#' @return graph with new edge weights
#' @export
#' @examples
#' graph.wt <- as(matrix(c(0,0,0,1,0,0,0,1,0), 3), "graphNEL")
#' graph.mt <- resample_edge_weights(graph.wt)
resample_edge_weights <- function(g, lB = -1, uB = 1, tp = 1) {
    if (length(uB) == 1) { uB <- c(0, uB) }
    if (length(lB) == 1) { lB <- c(lB, 0) }
    if (lB[1] - lB[2] == 0) {
        lB <- uB
    }
    if (uB[1] - uB[2] == 0) {
        uB <- lB
    }
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
#' @importFrom methods as
#' @importFrom MASS glm.nb
fulllin <- function(g1, d1, g2, d2, conf = TRUE, diff = 1,
                    errDist = "normal", ...) {
    mat1 <- as(g1, "matrix")
    mat2 <- as(g2, "matrix")
    mat1[which(mat1 != 0)] <- 1
    mat2[which(mat2 != 0)] <- 1
    dagtc <- nem::transitive.closure(mat1, mat=TRUE)
    df <- rbind(d1, d2)
    colnames(df) <- paste0("X", seq_len(ncol(df)))
    df <- as.data.frame(cbind(df, N = c(rep(1, nrow(d1)), rep(0, nrow(d2)))))

    n <- length(nodes(g1))

    if (diff) {
        dce <- mat1*0
        for (i in seq_len(n)) {
            for (j in seq_len(n)) {
                if (dagtc[i, j] == 1 & i != j) {
                    Z <- which(mat1[, i] == 1)
                    NX <- df[, i] * df$N
                    NZ <- df[, Z] * df$N
                    X <- df[, i]
                    Y <- df[, j]
                    Z <- df[, Z]
                    if (errDist %in% "normal") {
                        if (length(Z) > 0 & conf) {
                            C <- cov(cbind(Y, NX, NZ, X, Z))
                        } else {
                            C <- cov(cbind(Y, NX, X))
                        }
                        if (Matrix::rankMatrix(C[2:nrow(C), 2:ncol(C)]) <
                            nrow(C[2:nrow(C), 2:ncol(C)])) {
                            betas <- Gsolve(
                                C[2:nrow(C), 2:ncol(C)], C[2:nrow(C), 1]
                            )
                        } else {
                            betas <- solve(
                                C[2:nrow(C), 2:ncol(C)], C[2:nrow(C), 1]
                            )
                        }
                        dce[i, j] <- betas[1]
                    } else if (errDist %in% "nbinom") {
                        if (length(Z) > 0 & conf) {
                            NZ <- as(NZ, "matrix")
                            Z <- as(Z, "matrix")
                            betas <- glm.nb(Y ~ NX + df$N + X + NZ + Z, link = "identity")
                        } else {
                            betas <- glm.nb(Y ~ NX + df$N + X, link = "identity")
                        }
                        dce[i, j] <- betas$coefficients[2]
                    }
                }
            }
        }
    }
    res <- list(dce = dce, graph = g1)
    class(res) <- "dce"
    return(res)
}
