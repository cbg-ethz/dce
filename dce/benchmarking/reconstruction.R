library(tidyverse)

source("../R/utils.R")


# create graph
# df.graph <- data.frame(
#   from=c("A"),
#   to=c("B")
# )
# g <- tidygraph::as_tbl_graph(df.graph)
# ig <- igraph::as.igraph(g)
# graph <- igraph::igraph.to.graphNEL(ig)
# w <- graph@edgeData@data
# for (i in 1:length(w)) {
#   w[[i]]$weight <- 2.5 #runif(1, -2, 2)
# }
# graph@edgeData@data <- w

graph <- create_random_DAG(30, .2, c(-1,-.1), c(.1, 1))

adj <- t(as(graph, "matrix"))


# is graph topologically ordered?
nonZeros <- which(adj != 0, arr.ind=TRUE)

if (nrow(nonZeros) > 0) {
  if (any(nonZeros[, 1] - nonZeros[, 2] < 0) || any(diag(adj) != 0))
    stop("Input DAG must be topologically ordered!")
}


# start benchmark
set.seed(42)

sample.num <- 500
repetition.num <- 50
noise.list <- c(0.001, 0.05, 0.1, 0.5, 1, 2, 50)

df.bench <- purrr::map_df(rep(noise.list, each=repetition.num), function(noise.sd) {
  # generate data
  p <- length(nodes(graph))

  X <- matrix(0, nrow=sample.num, ncol=p)
  colnames(X) <- nodes(graph)
  X[, 1] <- 1 + rnorm(sample.num, mean=0, sd=noise.sd) # activate root node

  for (j in 2:p) {
    ij <- 1:(j - 1)
    X[, j] <- X[, j] + X[, ij, drop=FALSE] %*% adj[j, ij] + rnorm(sample.num, mean=0, sd=noise.sd)
  }


  # compute effects
  cor.mat <- cor(X)
  cov.mat <- cov(X)

  cau.mat <- purrr::map_dfc(1:p, function(x) {
    pcalg::idaFast(x, 1:p, cov.mat, graph)
  })


  # measure accuracy
  rmse <- function(m, o) {
    sqrt(mean(as.matrix((m - o)^2)))
  }
  cor.p <- function(m, o) {
    m.tmp <- as.matrix(m)
    o.tmp <- as.matrix(o)

    # diag(m.tmp) <- 0
    # diag(o.tmp) <- 0
    #
    # # flatten matrices
    # dim(m.tmp) <- NULL
    # dim(o.tmp) <- NULL

    nonzero.idx <- which(adj != 0, arr.ind=TRUE)

    cor(m.tmp[nonzero.idx], o.tmp[nonzero.idx], method="spearman")
  }

  data.frame(
    method=c(
      "correlation",
      "causal_effect"
    ),
    performance=c(
      "cor.p",
      "cor.p"
    ),
    value=c(
      cor.p(adj, cor.mat),
      cor.p(adj, cau.mat)
    )
  ) %>%
    mutate(noise.sd=as.factor(noise.sd))
})

df.bench %>%
  head

# plot
p.graph <- GGally::ggnet2(
  t(adj),
  label=TRUE,# edge.label="weights",
  arrow.size=12, arrow.gap=0.025
)

p.bench <- ggplot(df.bench, aes(x=noise.sd, y=value, color=method)) +
  geom_boxplot() +
  facet_grid( ~ performance) +
  theme_minimal()

cowplot::plot_grid(p.graph, p.bench, rel_widths=c(1,2))
ggsave("benchmarking.pdf")
