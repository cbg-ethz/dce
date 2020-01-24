#' @export
simulate_data <- function(
  dag, n = 100,
  dist.mean = 1000, dist.dispersion = 100
) {
  p <- length(nodes(dag))
  adj.mat <- igraph::as_adjacency_matrix(
    igraph::igraph.from.graphNEL(dag),
    attr="weight"
  ) %>%
    as.matrix %>%
    t

  # sanity checks
  stopifnot(is(dag, "graph"), p >= 2)

  if (any(adj.mat < 0)) {
    stop("Negative edge weights are not allowed!")
  }

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
  for (j in seq(2, p)) {
    ij <- seq_len(j - 1)
    betas <- adj.mat[j, ij]

    if (any(betas != 0)) {
      # current node has parents
      mu <- X[, ij, drop = FALSE] %*% betas
      X[, j] <- rnbinom(n, size=dist.dispersion, mu=mu)
    }
  }
  return(X)
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
simulate_data_old <- function (
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
                betas <- weightMatrix[j, ij]
                betas[betas < 0] <- 0
                mu <- X[, ij, drop = FALSE] %*% betas
                if (any(mu == 0)) {
                    mu <- rep(normpars[1], length(mu))
                }
                X[, j] <- rnegbin(n, mu = mu, theta = normpars[2])
            } else {
                X[, j] <- X[, j] + X[, ij, drop = FALSE] %*% weightMatrix[j,
                    ij]
            }
        }
        return(X)
    }
    else errMat
}
