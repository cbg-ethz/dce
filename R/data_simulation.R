#' Simulate data
#'
#' Generate data for given DAG.
#' @param graph Graph to simulate on
#' @param n Number of samples
#' @param dist_mean distribution mean as numeric
#' @param dist_dispersion distribution dispersion
#' (actually dispersion^-1) as a scalar
#' @param link special link function for the negative binomial
#' distribution
#' @param pop_size numeric for the population size, e.g., pop_size=1000 adds
#' 1000-n random genes not in the graph
#' @param latent number of latent variables
#' @param latent.fun uniform "unif" or exponential "exp" distribution of latent
#' coefficients
#' @return graph
#' @export
#' @rdname simulate_data-methods
#' @importFrom methods is
#' @importFrom pcalg wgtMatrix
#' @importFrom MASS rnegbin
#' @examples
#' dag <- create_random_DAG(30, 0.2)
#' X <- simulate_data(dag)
setGeneric(
    "simulate_data",
    function(
        graph, n = 100,
        dist_mean = 1000, dist_dispersion = 100,
        link = negative.binomial.special()$linkfun,
        pop_size = 0, latent = 0, latent.fun = "unif"
    ) {
        standardGeneric("simulate_data")
    },
    package = "dce"
)


setOldClass("igraph")
#' @rdname simulate_data-methods
#' @importFrom igraph V
#' @importFrom naturalsort naturalorder
setMethod(
    "simulate_data",
    signature = signature(graph = "igraph"),
    function(
        graph, n = 100,
        dist_mean = 1000, dist_dispersion = 100,
        link = negative.binomial.special()$linkfun,
        pop_size = 0, latent = 0, latent.fun = "unif"
    ) {
        mat <- as(igraph::as_adjacency_matrix(graph, attr = "weight"),
                  "matrix")
        colnames(mat) <- rownames(mat) <- V(graph)
        mat <- mat[naturalorder(rownames(mat)), naturalorder(colnames(mat))]
        simulate_data(mat,
            n, dist_mean, dist_dispersion,
            link, pop_size, latent
        )
    }
)


#' @rdname simulate_data-methods
setMethod(
    "simulate_data",
    signature = signature(graph = "graphNEL"),
    function(
        graph, n = 100,
        dist_mean = 1000, dist_dispersion = 100,
        link = negative.binomial.special()$linkfun,
        pop_size = 0, latent = 0, latent.fun = "unif"
    ) {
        a <- as(graph, "matrix")
        a <- a[naturalorder(rownames(a)), naturalorder(colnames(a))]
        simulate_data(
            a,
            n, dist_mean, dist_dispersion,
            link, pop_size, latent
        )
    }
)


#' @rdname simulate_data-methods
#' @importFrom naturalsort naturalorder
setMethod(
    "simulate_data",
    signature = signature(graph = "matrix"),
    function(
        graph, n = 100,
        dist_mean = 1000, dist_dispersion = 100,
        link = negative.binomial.special()$linkfun,
        pop_size = 0, latent = 0, latent.fun = "unif"
    ) {
        start <- 2
        p <- dim(graph)[[1]]
        if (latent > 0) {
            if (latent.fun == "unif") {
                H1 <- matrix(runif(p * latent, -1, 1), latent, p)
            } else if (latent.fun == "exp") {
                H1 <- matrix(exp(p * latent), latent, p)
            }
            H0 <- matrix(0, p + latent, latent)
            graph <- cbind(H0, rbind(H1, graph))
            start <- latent + 1
            p <- dim(graph)[[1]]
        }
        if (pop_size < p) {
            pop_size <- p
        }

        # sanity checks
        stopifnot(p >= 2)

        nonzero_idx <- which(graph != 0, arr.ind = TRUE)
        if (
            (nrow(nonzero_idx) > 0) &&
            (
                any(nonzero_idx[, 2] - nonzero_idx[, 1] < 0) ||
                any(diag(graph) != 0)
            )
        ) {
            stop("Input DAG must be topologically ordered!")
        }

        # setup data
        X <- matrix(rnbinom(n * p,
                            size = dist_dispersion,
                            mu = dist_mean),
                    nrow = n, ncol = p)
        colnames(X) <- colnames(graph)

        # simulate data
        for (j in seq(start, p)) {
            ij <- seq_len(j - 1)
            betas <- graph[ij, j]

            if (any(betas != 0)) {
                # current node has parents
                mu <- link(X[, ij, drop = FALSE] %*% betas, offset = 1)
                X[, j] <- rpois(n, lambda = mu)
            }
        }

        X <- X[, naturalorder(colnames(X))]
        if (pop_size > p & latent == 0) {
            Y <- matrix(rnbinom(n * (pop_size - p),
                                size = dist_dispersion,
                                mu = dist_mean),
                        nrow = n, ncol = pop_size - p)
            colnames(Y) <- paste0("n", (p + 1):pop_size)
            X <- cbind(X, Y)
        } else if (pop_size > p & latent > 0) {
            H <- matrix(
                runif(latent * (pop_size - p + latent), -1, 1),
                latent, pop_size - p + latent
            )
            Y <- X[, seq_len(latent), drop = FALSE] %*% H
            Y <- apply(Y, 2, function(x) {
                y <- rpois(n, lambda = link(x, offset = 1))
                return(y)
            })
            colnames(Y) <- paste0("n", (p + 1 - latent):pop_size)
            X <- cbind(X, Y)
        }
        if (latent > 0) {
            X <- X[, -seq_len(latent), drop = FALSE]
        }
        return(X)
    }
)
