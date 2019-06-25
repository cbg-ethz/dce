library(tidyverse)
library(purrr)

library(graph)
library(pcalg)
library(assertthat)

compute_causal_effects <- function(graph, df.expr) {
  node.list <- nodes(graph)

  cov.mat <- cov(df.expr)
  assertthat::are_equal(dim(cov.mat)[[1]], length(node.list))

  n <- length(node.list)
  purrr::map_dfc(1:n, function(x) {
    pcalg::idaFast(x, 1:n, cov.mat, graph)
  })
}

compute_differential_causal_effects <- function(
  graph.ctrl, df.expr.ctrl,
  graph.mut, df.expr.mut
) {
  ce.ctrl <- compute_causal_effects(graph.ctrl, df.expr.ctrl)
  ce.mut <- compute_causal_effects(graph.mut, df.expr.mut)

  as.matrix(ce.ctrl - ce.mut)
}
