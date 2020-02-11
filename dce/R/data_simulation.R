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
simulate_data <- function(
  dag, n = 100,
  dist.mean = 1000, dist.dispersion = 100,
  link.log.base=exp(1)
) {
  p <- length(nodes(dag))
  adj.mat <- igraph::as_adjacency_matrix(
    igraph::igraph.from.graphNEL(dag),
    attr=if (numEdges(dag) > 0) "weight" else NULL
  ) %>%
    as.matrix %>%
    t

  # sanity checks
  stopifnot(is(dag, "graph"), p >= 2)

  nonzero.idx <- which(adj.mat != 0, arr.ind = TRUE)
  if (nrow(nonzero.idx) > 0) {
    if (any(nonzero.idx[, 1] - nonzero.idx[, 2] < 0) || any(diag(adj.mat) != 0)) {
      stop("Input DAG must be topologically ordered!")
    }
  }

  # setup data
  X <- matrix(rnbinom(n * p, size=dist.dispersion, mu=dist.mean), nrow=n, ncol=p)
  colnames(X) <- nodes(dag)

  # simulate data
  link <- make.log.link(link.log.base)

  for (j in seq(2, p)) {
    ij <- seq_len(j - 1)
    betas <- adj.mat[j, ij]

    if (any(betas != 0)) {
      # current node has parents
      if (link.log.base == 0) {
        betasX <- X[, ij, drop = FALSE] %*% betas
        mu <- dist.mean + betasX - min(betasX)
      } else {
        mu <- link$linkinv(log(dist.mean, link.log.base) + scale(X[, ij, drop = FALSE], scale=FALSE) %*% betas)
      }
      X[, j] <- rnbinom(n, size=dist.dispersion, mu=mu)
    }
  }
  return(X)
}
