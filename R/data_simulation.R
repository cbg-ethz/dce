#' Simulate data
#'
#' Generate data for given DAG. The flexible framework allows for different
#' distributions for source and child nodes. Default distributions are
#' negative binomial (with mean = 100 and 1/dispersion = 100), and poisson,
#' respectively.
#' @param graph Graph to simulate on
#' @param n Number of samples
#' @param dist_fun distribution function for nodes without parents
#' @param dist_args list of arguments for dist_fun
#' @param child_fun distribution function for nodes with parents
#' @param child_args list of arguments for child_fun
#' @param child_dep link_fun computes an output for the expression of
#' nodes without parents. this output is than used as input for child_fun.
#' child_dep defines the parameter (a a string) of child_fun, which is used for
#' the input. E.g., the link_fun is the identity and the child_fun is rnorm, we
#' usually set child_dep = "mean".
#' @param link_fun special link function for the negative binomial
#' distribution
#' @param link_args list of arguments for link_fun
#' @param pop_size numeric for the population size, e.g., pop_size=1000 adds
#' 1000-n random genes not in the graph
#' @param latent number of latent variables
#' @param latent_fun uniform "unif" or exponential "exp" distribution of latent
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
        dist_fun = rnbinom, dist_args = list(mu = 1000, size = 100),
        child_fun = rpois, child_args = list(), child_dep = "lambda",
        link_fun = negative.binomial.special()$linkfun,
        link_args = list(offset = 1),
        pop_size = 0, latent = 0, latent_fun = "unif"
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
        dist_fun = rnbinom, dist_args = list(mu = 1000, size = 100),
        child_fun = rpois, child_args = list(), child_dep = "lambda",
        link_fun = negative.binomial.special()$linkfun,
        link_args = list(offset = 1),
        pop_size = 0, latent = 0, latent_fun = "unif"
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
        dist_fun = rnbinom, dist_args = list(mu = 1000, size = 100),
        child_fun = rpois, child_args = list(), child_dep = "lambda",
        link_fun = negative.binomial.special()$linkfun,
        link_args = list(offset = 1),
        pop_size = 0, latent = 0, latent_fun = "unif"
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
        dist_fun = rnbinom, dist_args = list(mu = 1000, size = 100),
        child_fun = rpois, child_args = list(), child_dep = "lambda",
        link_fun = negative.binomial.special()$linkfun,
        link_args = list(offset = 1),
        pop_size = 0, latent = 0, latent_fun = "unif"
    ) {
        start <- 2
        p <- dim(graph)[[1]]
        if (latent > 0) {
            if (latent_fun == "unif") {
                H1 <- matrix(runif(p * latent, -1, 1), latent, p)
            } else if (latent_fun == "exp") {
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
        X <- matrix(do.call(dist_fun, c(n * p, dist_args)),
                    nrow = n, ncol = p)
        colnames(X) <- colnames(graph)

        # simulate data
        for (j in seq(start, p)) {
            ij <- seq_len(j - 1)
            betas <- graph[ij, j]

            if (any(betas != 0)) {
                # current node has parents
                mu <- do.call(link_fun, c(list(X[, ij, drop = FALSE] %*% betas),
                                          link_args))
                X[, j] <- unlist(lapply(mu, function(x) {
                    y_list <- c(1, child_args)
                    y_list[[child_dep]] <- x
                    y <- do.call(child_fun, y_list)
                    return(y)
                }))
            }
        }

        X <- X[, naturalorder(colnames(X))]
        if (pop_size > p & latent == 0) {
            Y <- matrix(do.call(dist_fun, c(n * (pop_size - p), dist_args)),
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
                y <- do.call(child_fun,
                             c(list(n, do.call(link_fun, c(x, link_args)))))
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
