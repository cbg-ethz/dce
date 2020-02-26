#' Simulate data
#'
#' Generate data for given DAG.
#' @param dag Graph to simulate on
#' @param n Number of samples
#' @param dist.mean distribution mean as numeric
#' @param dist.dispersion distribution dispersion
#' (actually dispersion^-1) as a scalar
#' @return graph
#' @export
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
        dist.mean = 1000, dist.dispersion = 100,
        link = negative.binomial.special()$linkfun
    ) {
        standardGeneric("simulate_data")
    },
    package = "dce"
)


setOldClass("igraph")
setMethod(
    "simulate_data",
    signature = signature(graph = "igraph"),
    function(
        graph, n = 100,
        dist.mean = 1000, dist.dispersion = 100,
        link = negative.binomial.special()$linkfun
    ) {
        simulate_data(
            as(igraph::as_adjacency_matrix(
                graph,
                attr = if ("weight" %in% igraph::edge_attr_names(graph)) "weight" else NULL
            ), "matrix"),
            n, dist.mean, dist.dispersion,
            link
        )
    }
)


setMethod(
    "simulate_data",
    signature = signature(graph = "graphNEL"),
    function(
        graph, n = 100,
        dist.mean = 1000, dist.dispersion = 100,
        link = negative.binomial.special()$linkfun
    ) {
        simulate_data(
            as(graph, "matrix"),
            n, dist.mean, dist.dispersion,
            link
        )
    }
)


setMethod(
    "simulate_data",
    signature = signature(graph = "matrix"),
    function(
        graph, n = 100,
        dist.mean = 1000, dist.dispersion = 100,
        link = negative.binomial.special()$linkfun
    ) {
        p <- dim(graph)[[1]]

        # sanity checks
        stopifnot(p >= 2)

        nonzero.idx <- which(graph != 0, arr.ind = TRUE)
        if (
            (nrow(nonzero.idx) > 0) &&
            (
                any(nonzero.idx[, 2] - nonzero.idx[, 1] < 0) ||
                any(diag(graph) != 0)
            )
        ) {
            stop("Input DAG must be topologically ordered!")
        }

        # setup data
        X <- matrix(rnbinom(n * p, size = dist.dispersion, mu = dist.mean), nrow = n, ncol = p)
        colnames(X) <- colnames(graph)

        # simulate data
        for (j in seq(2, p)) {
          ij <- seq_len(j - 1)
          betas <- graph[ij, j]

          if (any(betas != 0)) {
            # current node has parents
            mu <- link(X[, ij, drop = FALSE] %*% betas, offset = dist.mean*0 + 1)
            X[, j] <- rnbinom(n, size = dist.dispersion, mu = mu)
          }
        }
        return(X)
    }
)
