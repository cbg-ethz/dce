#' Graph to DAG
#'
#' Converts a general graph to a dag with minimum distance to the
#' original graph. The general idea is to transitively close the graph
#' to detect cycles and remove them based on the rule "the more outgoing
#' edges a node has, the more likely it is that incoming edges from a
#' cycle will be deleted, and vice versa. However, this is too rigorous
#' and deletes too many edges, which do not lead to a cycle. These
#' edges are added back in the final step.
#' @param g graph as adjacency matrix
#' @param tc if TRUE computes the transitive closure
#' @author Ken Adams
#' @return dag as adjacency matrix
#' @export
#' @examples
#' g <- matrix(c(1,0,1,0, 1,1,0,0, 0,1,1,0, 1,1,0,1), 4, 4)
#' rownames(g) <- colnames(g) <- LETTERS[1:4]
#' dag <- g2dag(g)
g2dag <- function(g, tc = FALSE) {
    ord <- order(apply(g, 1, sum) - apply(g, 2, sum), decreasing = 1)
    g <- g[ord, ord]
    cyc <- intersect(which(g+t(g) > 1), which(lower.tri(g) == TRUE))
    g[cyc] <- 0
    diag(g) <- 1
    g2 <- g
    g0 <- g*0
    epiNEM::HeatmapOP(g, Colv = 0, Rowv = 0)
    for (i in seq_len(nrow(g))) {
        g2 <- g2%*%g
        g2[which(g2 > 0)] <- 1
        ord <- order(apply(g2, 1, sum) - apply(g2, 2, sum), decreasing = 1)
        g2 <- g2[ord, ord]
        g <- g[ord, ord]
        cyc <- intersect(which(g2+t(g2) > 1), which(lower.tri(g2) == TRUE))
        g2[cyc] <- 0
        if (all(g0 == g2)) { break() }
        g0 <- g2
        epiNEM::HeatmapOP(g2, Colv = 0, Rowv = 0)
    }
    if (!tc) {
        print(g)
        g3 <- g2*0
        g3[which(g2 == 1 & g == 1)] <- 1
        for (e in which(g == 1 & g3 == 0)) {
            for (a in 1:ncol(g)) {
                b <- e - a*nrow(g)
                if (b < 0) {
                    b <- e - (a-1)*nrow(g)
                    break()
                }
            }
            print(e)
            print(b)
            print(a)
            if (sum(g[b, ]) >= sum(g[a, ])) {
                g3[b, a] <- 1
            }
            g4 <- mytc(g3)
            ord <- order(apply(g4, 1, sum) - apply(g4, 2, sum), decreasing = 1)
            g4 <- g4[ord, ord]
            if (any(g2+t(g2) > 1)) {
                g3[a, b] <- 0
            }
        }
        g2 <- g3
    }
    return(g2)
}
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
    node.labels = paste0("n", as.character(seq_len(n)))
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
#' @importFrom zetadiv glm.cons
fulllin <- function(g1, d1, g2, d2, conf = TRUE,
                    errDist = "nbinom", theta = NULL,
                    link.log.base = exp(1),
                    partial = TRUE, ...) {
    mat1 <- as(g1, "matrix")
    mat2 <- as(g2, "matrix")
    mat1[which(mat1 != 0)] <- 1
    mat2[which(mat2 != 0)] <- 1
    dagtc <- nem::transitive.closure(mat1, mat=TRUE)
    df <- rbind(d1, d2)
    if (is.null(theta)) {
        theta <- estimateTheta(df)
    }
    colnames(df) <- paste0("X", seq_len(ncol(df)))
    df <- as.data.frame(cbind(df, N = c(rep(1, nrow(d1)), rep(0, nrow(d2)))))
    n <- length(nodes(g1))
    dce <- mat1*0
    dce.p <- mat1*NA
    glmfun <- function(formula) {
        fun <- "glm2"
        
        if (link.log.base == 0) {
            link <- "identity"
        } else {
            link <- make.log.link(link.log.base)
        }

        if (fun %in% "glm.nb") {
            fit <- MASS::glm.nb(formula, link = link, ...)
        } else if (fun %in% "glm2") {
            fit <- glm2::glm2(formula, family = MASS::negative.binomial(
                                                theta=theta,
                                                link=link), ...)
        } else if (fun %in% "glm.cons") {
            fit <- zetadiv::glm.cons(formula,
                                     family = MASS::negative.binomial(
                                                        theta=theta,
                                                        link=link),
                                     cons = 1, ...)
        } else if (fun %in% "gauss") {
            fit <- glm(formula, family = "gaussian", ...)
        }
        return(fit)
    }
    if (partial) {
        for (i in seq_len(n)) {
            Xidx <- which(mat1[, i] == 1)
            if (length(Xidx) > 0) {
                Y <- df[, i]
                NX <- as(df[, Xidx] * df$N, "matrix")
                X <- as(df[, Xidx], "matrix")
                N <- df$N
                if (errDist %in% "normal") {
                    C <- cov(cbind(Y, NX, N, X))
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
                    dce[Xidx, i] <- betas[1:length(Xidx)]
                } else if (errDist %in% "nbinom") {
                    fit <- glmfun(Y ~ NX + N + X)
                    coef.mat <- summary(fit)$coefficients
                    dce[Xidx, i] <- coef.mat[2:(length(Xidx)+1), 1]
                    dce.p[Xidx, i] <- coef.mat[2:(length(Xidx)+1), 4]
                }
            }
        }
    } else {
        for (i in seq_len(n)) {
            for (j in seq_len(n)) {
                if (dagtc[i, j] == 1 & i != j) {
                    Z <- which(mat1[, i] == 1)
                    NX <- df[, i] * df$N
                    NZ <- as(df[, Z] * df$N, "matrix")
                    N <- df$N
                    X <- df[, i]
                    Y <- df[, j]
                    Z <- as(df[, Z], "matrix")
                    if (errDist %in% "normal") {
                        if (length(Z) > 0 & conf) {
                            C <- cov(cbind(Y, NX, NZ, N, X, Z))
                        } else {
                            C <- cov(cbind(Y, NX, N, X))
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
                            fit <- glmfun(Y ~ NX + N + X + NZ + Z)
                        } else {
                            fit <- glmfun(Y ~ NX + N + X)
                        }
                        coef.mat <- summary(fit)$coefficients
                        dce[i, j] <- coef.mat[2, 1]
                        dce.p[i, j] <- coef.mat[2, 4]
                    }
                }
            }
        }
    }
    res <- list(dce = dce, dce.p = dce.p, graph = g1)
    class(res) <- "dce"
    return(res)
}
#' @importFrom edgeR DGEList calcNormFactors estimateDisp
#' @noRd
estimateTheta <- function(data) {
    if (ncol(data) < 5) {
        mus <- apply(data, 2, mean)
        sigmas <- apply(data, 2, sd)
        thetas <- mus^2/(sigmas^2 - mus)
        theta <- mean(thetas)
    } else {
        y <- edgeR::DGEList(counts=t(data))
        y <- edgeR::calcNormFactors(y)
        y <- edgeR::estimateDisp(y)
        theta <- 1/y$common.dispersion
    }
    return(theta)
}

#' @export
make.log.link <- function(base=exp(1)) {
    structure(list(
        linkfun=function(mu) { log(mu, base) },
        linkinv=function(eta) { base**(eta) },
        mu.eta=function(eta) { base**(eta) },
        valideta=function(eta) { TRUE }
    ), class="link-glm")
}
