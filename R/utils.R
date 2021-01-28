#' Permutation test for (partial) correlation on non-Gaussian data
#'
#' Computes the significance of (partial) correlation based
#' on permutations of the observations
#' @param x wild type data set
#' @param y mutant data set
#' @param iter number of iterations (permutations)
#' @param fun function to compute the statistic, e.g., cor or pcor
#' @param ... additional arguments for function 'fun'
#' @return matrix of p-values
#' @export
#' @examples
#' x <- matrix(rnorm(100),10,10)
#' y <- matrix(rnorm(100),10,10)
#' permutation_test(x,y,iter=10)
permutation_test <- function(x, y, iter = 1000, fun = pcor, ...) {
    z <- fun(y, ...) - fun(x, ...)
    p <- z * 0
    for (i in seq_len(iter)) {
        xy <- rbind(x, y)
        xyp <- xy[sample(seq_len(nrow(xy)), nrow(xy)), ]
        xp <- xyp[seq_len(nrow(x)), ]
        yp <- xyp[-seq_len(nrow(x)), ]
        zp <- fun(yp, ...) - fun(xp, ...)
        idx <- which(abs(zp) >= abs(z))
        p[idx] <- p[idx] + 1
    }
    p <- p / iter
    return(p)
}

#' graph to adjacency
#'
#' From graphNEL with 0 edge weights to proper adjacency matrix
#' @param g graphNEL object
#' @export
#' @return graph as adjacency matrix
#' @examples
#' dag <- create_random_DAG(30, 0.2)
#' adj <- as_adjmat(dag)
as_adjmat <- function(g) {
    adj <- as(g, "matrix")
    for (p in names(g@edgeData@data)) {
        a <- gsub("\\|.*", "", p)
        b <- gsub(".*\\|", "", p)
        adj[a, b] <- 1
    }
    return(adj)
}

#' Partial correlation
#'
#' Robust partial correlation of column variables of a
#' numeric matrix
#' @param x matrix
#' @param ... additional arguments for functions 'cor'
#' @importFrom ppcor pcor
#' @export
#' @return matrix of partial correlations
#' @examples
#' x <- matrix(rnorm(100),10,10)
#' pcor(x)
pcor <- function(x, ...) {
    rho <- try(ppcor::pcor(x, ...), silent = TRUE)
    if (length(grep("Error", rho)) > 0) {
        warning(paste0("Moore-Penrose generalized matrix invers in ",
                       "function ppcor::pcor crashed. Using ",
                       "MAAS::ginv instead."))
        omega <- cor(x, ...)
        p <- Gsolve(omega)
        pdiag <- diag(p) %*% t(diag(p))
        rho <- -p / (pdiag^0.5)
        diag(rho) <- 1
        rho[which(is.na(rho) | is.infinite(rho))] <- 0
    } else {
        rho <- rho$estimate
    }

    rownames(rho) <- colnames(rho) <- colnames(x)

    return(rho)
}


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
#' @importFrom mnem transitive.closure
#' @examples
#' g <- matrix(c(1,0,1,0, 1,1,0,0, 0,1,1,0, 1,1,0,1), 4, 4)
#' rownames(g) <- colnames(g) <- LETTERS[seq_len(4)]
#' dag <- g2dag(g)
g2dag <- function(g, tc = FALSE) {
    ord <- order(apply(g, 1, sum) - apply(g, 2, sum), decreasing = 1)
    g <- g[ord, ord]
    cyc <- intersect(which(g + t(g) > 1), which(lower.tri(g) == TRUE))
    g[cyc] <- 0
    diag(g) <- 1
    g2 <- g
    g0 <- g * 0
    epiNEM::HeatmapOP(g, Colv = 0, Rowv = 0)
    for (i in seq_len(nrow(g))) {
        g2 <- g2 %*% g
        g2[which(g2 > 0)] <- 1
        ord <- order(apply(g2, 1, sum) - apply(g2, 2, sum), decreasing = 1)
        g2 <- g2[ord, ord]
        g <- g[ord, ord]
        cyc <- intersect(which(g2 + t(g2) > 1), which(lower.tri(g2) == TRUE))
        g2[cyc] <- 0
        if (all(g0 == g2)) {
            break
        }
        g0 <- g2
        epiNEM::HeatmapOP(g2, Colv = 0, Rowv = 0)
    }
    if (!tc) {
        print(g)
        g3 <- g2 * 0
        g3[which(g2 == 1 & g == 1)] <- 1
        for (e in which(g == 1 & g3 == 0)) {
            for (a in seq_len(ncol(g))) {
                b <- e - a * nrow(g)
                if (b < 0) {
                    b <- e - (a - 1) * nrow(g)
                    break
                }
            }
            print(e)
            print(b)
            print(a)
            if (sum(g[b, ]) >= sum(g[a, ])) {
                g3[b, a] <- 1
            }
            g4 <- transitive.closure(g3)
            ord <- order(apply(g4, 1, sum) - apply(g4, 2, sum), decreasing = 1)
            g4 <- g4[ord, ord]
            if (any(g2 + t(g2) > 1)) {
                g3[a, b] <- 0
            }
        }
        g2 <- g3
    }
    return(g2)
}


#' @importFrom MASS ginv
#' @noRd
Gsolve <- function(a, b=NULL, ...) {
    if (!is.null(b)) {
        x <- t(solve(diag(1, nrow(a)), MASS::ginv(a, ...) %*% b))
    } else {
        x <- MASS::ginv(a, ...)
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
        sign(truth) == sign(inferred) & inferred != 0,
        na.rm = TRUE
    )
    fn <- sum(
        abs(truth) > cutoff &
        abs(inferred) <= cutoff &
        inferred != 0,
        na.rm = TRUE
    )
    fp <- sum(
        abs(truth) <= cutoff &
        abs(inferred) > cutoff |
        (
            abs(truth) > cutoff &
            abs(inferred) > cutoff &
            sign(truth) != sign(inferred)
        ) &
        inferred != 0,
        na.rm = TRUE
    )
    tn <- sum(
        abs(truth) <= cutoff &
        abs(inferred) <= cutoff &
        inferred != 0,
        na.rm = TRUE
    )

    return(data.frame(
        true.positive = tp,
        false.positive = fp,
        true.negative = tn,
        false.negative = fn
    ))
}


#' Create random DAG (topologically ordered)
#'
#' Creates a DAG according to given parameters.
#' @param node_num Number of nodes
#' @param prob Probability of creating an edge
#' @param eff_min Lower bound for edge weights
#' @param eff_max Upper bound for edge weights
#' @param node_labels Node labels
#' @author Martin Pirkl
#' @return graph
#' @export
#' @importFrom methods new
#' @examples
#' dag <- create_random_DAG(30, 0.2)
create_random_DAG <- function(
    node_num, prob,
    eff_min = -1, eff_max = 1,
    node_labels = paste0("n", as.character(seq_len(node_num)))
) {
    stopifnot(
        node_num >= 2,
        is.numeric(prob),
        length(prob) == 1,
        0 <= prob,
        prob <= 1,
        is.numeric(eff_min),
        is.numeric(eff_max),
        eff_min <= eff_max
    )

    # create (directed) adjacency matrix
    mat <- matrix(rbinom(node_num * node_num, 1, prob), node_num, node_num)
    mat[lower.tri(mat)] <- 0
    diag(mat) <- 0

    # assign effects
    mat[mat != 0] <- runif(sum(mat != 0), min = eff_min, max = eff_max)

    # set node labels
    rownames(mat) <- colnames(mat) <- node_labels

    # topologically sort graph
    nodes_sorted <- igraph::topo_sort(
        igraph::graph_from_adjacency_matrix(mat, weighted = TRUE)
    )
    mat <- mat[nodes_sorted, nodes_sorted]

    # return graph
    igraph::as_graphnel(
        igraph::graph_from_adjacency_matrix(mat, weighted = TRUE)
    )
}


#' Resample network edge weights
#'
#' Takes a graph and modifies edge weights.
#' @param g original graph
#' @param tp fraction of edge weights which will be modified
#' @param mineff minimal differential effect size
#' @param maxeff maximum effect effect size or standard deviation,
#' if method is "gauss"
#' @param method method for drawing the differential for the causal effects.
#' Can be "runif", "exp" or "gauss".
#' @author Martin Pirkl
#' @return graph with new edge weights
#' @export
#' @examples
#' graph.wt <- as(matrix(c(0,0,0,1,0,0,0,1,0), 3), "graphNEL")
#' graph.mt <- resample_edge_weights(graph.wt)
resample_edge_weights <- function(g, tp = 0.5,
                                  mineff = 1, maxeff = 2,
                                  method = "runif") {
    gold <- g
    g <- as(g, "matrix")
    changes <- floor((1 - tp) * sum(g != 0))
    gbin <- g
    gbin[which(gbin != 0)] <- 1
    gtr <- gtr2 <- mnem::transitive.reduction(gbin)
    diag(gtr) <- 1
    noedges <- 1
    start <- TRUE
    keep <- NULL
    while (changes > 0 & sum(gtr2) != 0 | start) {
        start <- FALSE
        if (changes == 1 & sum(gtr2) == 1) {
            keep <- c(keep, which(gtr2 == 1))
        } else {
            keep <- c(keep,
                      sample(which(gtr2 == 1),
                             min(changes, sum(gtr2))))
        }
        changes <- changes - sum(gtr2)
        diag(gtr2) <- 1
        gtr2[noedges] <- 1
        gtr2[keep] <- 1
        gtr2 <- gtr2 %*% gtr
        gtr2[which(gtr2 > 0)] <- 1
        noedges <- which(gbin == 0 & gtr2 == 1)
        gtr2[noedges] <- 0
        gtr2[keep] <- 0
    }
    g2 <- g
    change <- which(g != 0)
    change <- unlist(lapply(change, function(x) {
        eff <- g2[x]
        die <- sample(c(0, 1), 1)
        if (die) {
            if (method == "runif") {
                effnew <- runif(1, eff + mineff, eff + maxeff)
            } else if (method == "gauss") {
                effnew <- eff + mineff + abs(rnorm(1, 0, maxeff))
            } else if (method == "exp") {
                effnew <- eff + mineff + rexp(1, maxeff)
            }
        } else {
            if (method == "runif") {
                effnew <- runif(1, eff - maxeff, eff - mineff)
            } else if (method == "gauss") {
                effnew <- eff - mineff - abs(rnorm(1, 0, maxeff))
            } else if (method == "exp") {
                effnew <- eff - mineff - rexp(1, maxeff)
            }
        }
        return(effnew)
    }))
    g[which(g != 0)] <- change
    g[keep] <- g2[keep]
    edges <- which(gbin == 1, arr.ind = TRUE)
    for (i in seq(nrow(edges))) {
        edge <- paste0(
            rownames(gbin)[edges[i, 1]], "|", rownames(gbin)[edges[i, 2]]
        )
        gold@edgeData@data[[edge]]$weight <- g[edges[i, 1], edges[i, 2]]
    }
    return(gold)
}


#' @importFrom edgeR DGEList calcNormFactors estimateDisp
#' @noRd
estimateTheta <- function(data, ...) {
    if (is.null(dim(data))) {
        data <- matrix(data, length(data))
    }
    if (ncol(data) < 5) {
        mus <- apply(data, 2, mean)
        sigmas <- apply(data, 2, sd)
        thetas <- mus^2 / (sigmas^2 - mus)
        theta <- mean(thetas)
    } else {
        y <- edgeR::DGEList(counts = t(data))
        y <- edgeR::calcNormFactors(y)
        y <- edgeR::estimateDisp(y, ...)
        theta <- 1 / y$common.dispersion
    }
    return(theta)
}


#' @noRd
make.log.link <- function(base=exp(1)) {
    structure(list(
        linkfun = function(mu) log(mu, base),
        linkinv = function(eta) base**eta,
        mu.eta = function(eta) base**eta,
        valideta = function(eta) TRUE
    ), class = "link-glm")
}


#' Compute the true casual effects of a simulated dag
#'
#' This function takes a DAG with edgeweights as input and computes
#' the causal effects of all nodes on all direct and indirect children in the
#' DAG. Alternatively see pcalg::causalEffect for pairwise computation.
#' @param g graphNEL object
#' @param partial if FALSE computes the total causal effects and not just
#' the partial edge effects
#' @author Martin Pirkl
#' @return matrix of causal effects
#' @export
#' @importFrom pcalg causalEffect
#' @import graph tidyverse
#' @importFrom expm %^%
#' @examples
#' graph.wt <- as(matrix(c(0,0,0,1,0,0,0,1,0), 3), "graphNEL")
#' trueEffects(graph.wt)
trueEffects <- function(g, partial = FALSE) {
    a <- as(g, "matrix")
    if (partial) {
        ae <- a
    } else {
        ae <- a
        for (i in 2:nrow(a)) {
            ae <- ae + a %^% i
        }
    }
    return(ae)
}

#' @noRd
permutation_thresholding <- function(X, quantile=0.99) {
    N <- 50
    r <- min(dim(X))
    X <- scale(X)

    evals <- matrix(0, nrow = N, ncol = r)
    for (i in seq_len(N)) {
        X_perm <- apply(X, 2, function(xx) sample(xx))
        evals[i, ] <- svd(scale(X_perm))$d[seq_len(r)]
    }
    thresholds <- apply(evals, 2, function(xx) quantile(xx, probs = quantile))

    # last which crosses the threshold
    return(max(which(c(TRUE, svd(X)$d > thresholds))) - 1)
}

#' Estimate number of latent confounders
#' Compute the true casual effects of a simulated dag
#'
#' This function takes a DAG with edgeweights as input and computes
#' the causal effects of all nodes on all direct and indirect children in the
#' DAG. Alternatively see pcalg::causalEffect for pairwise computation.
#' @param X1 data matrix corresponding to the first condition
#' @param X2 data matrix corresponding to the second condition
#' @param method a string indicating the method used for estimating the number
#' of latent variables
#' @author Domagoj Ćevid
#' @return estimated number of latent variables
#' @export
#' @importFrom stats prcomp
#' @examples
#' graph1 <- create_random_DAG(node_num = 100, prob = .1)
#' graph2 <- resample_edge_weights(graph1, tp=0.15)
#' X1 <- simulate_data(graph1, n=200, latent = 3)
#' X2 <- simulate_data(graph2, n=200, latent = 3)
#' estimate_latent_count(X1, X2)
estimate_latent_count <- function(X1, X2, method = "auto") {
    if (method == "auto") {
        # this function estimates the number of principal components by
        # comparing the singular values to what they would be if the
        # variables were independent which is estimated by permuting the
        # columns of the data matrix
        permutation_thresholding <- function(X, quantile=0.99) {
            N <- 50
            r <- min(dim(X))
            X <- scale(X)

            evals <- matrix(0, nrow = N, ncol = r)
            for (i in seq_len(N)) {
                X_perm <- apply(X, 2, function(xx) sample(xx))
                evals[i, ] <- svd(scale(X_perm))$d[seq_len(r)]
            }
            thresholds <- apply(evals,
                                2, function(xx) quantile(xx, probs = quantile))

            # last which crosses the threshold
            return(max(which(c(TRUE, svd(X)$d > thresholds))) - 1)
        }

        q1 <- permutation_thresholding(X1)
        q2 <- permutation_thresholding(X2)

        return(ceiling((q1 + q2) / 2))
    }

    if (method == "kim") {
        # This function looks at the knee point of a scree plot.

        X <- cbind(X1, X2)
        fit_pca <- prcomp(scale(X))

        scree <- fit_pca$sdev
        scree <- scree[seq_len(round(length(scree)/2))]
        values <- seq(length(scree))

        d1 <- diff(scree) / diff(values) # first derivative
        d2 <- diff(d1) / diff(values[-1]) # second derivative
        idx <- which.max(abs(d2))

        return(idx)
    }

    if (method=="cluster") {
        X <- cbind(X1, X2)
        fit_pca <- prcomp(scale(X))

        scree <- fit_pca$sdev
        clust <- kmeans(scree,centers=c(scree[1],scree[2],
                                        scree[round(length(scree)/2+1)],
                                        scree[length(scree)]))$cluster
        idx <- sum(clust==clust[1])
        return(idx)
    }
}
