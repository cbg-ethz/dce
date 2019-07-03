library(tidyverse)
library(purrr)

library(graph)
library(pcalg)
library(assertthat)

#' Compute the causal effect of all nodes on al other nodes of a given DAG
#'
#' This function takes a DAG and gene expression data as input and computes
#' the causal effect of all nodes on all other nodes in the DAG.
#' @param graph DAG as a graphNEL object
#' @param df.expr genes expression profiles with rows as observations and
#' columns as variables (nodes)
#' @author Kim Jablonski
#' @return vector of causal effects
#' @export
#' @importFrom purrr map_dfc
#' @importFrom pcalg idaFast
#' @importFrom assertthat are_equal
#' @import graph tidyverse
#' @examples
#' dag <- matrix(c(0,0,0,1,0,0,0,1,0), 3)
#' colnames(dag) <- rownames(dag) <- seq_len(3)
#' dag <- as(dag, "graphNEL")
#' d <- matrix(rnorm(100*3), 100)
#' colnames(d) <- seq_len(3)
#' compute_causal_effects(dag, d)
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

  res <- t(as.matrix(ce.ctrl - ce.mut))
  class(res) <- "dce"
  return(res)
}

plot.dce <- function(x, ...) {
}
