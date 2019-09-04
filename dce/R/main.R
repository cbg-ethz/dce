library(tidyverse)
library(purrr)

library(graph)
library(pcalg)
library(assertthat)

#' Compute the causal effects
#'
#' This function takes a DAG and gene expression data as input and computes
#' the causal effects of all nodes on all other nodes in the DAG.
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

#' Compute the differential causal effects
#'
#' This function takes a DAG and gene expression data of two distinct types
#' of populations (e.g. tumor versus normal) as input and computes
#' the differential causal effects of all nodes on all other nodes in the DAG
#' between the tow populations of samples.
#' @param graph.ctrl DAG as a graphNEL object of the control population
#' @param df.expr.ctrl genes expression profiles with rows as observations and
#' columns as variables (nodes) of the control population
#' @param graph.mut DAG as a graphNEL object of the tumor population
#' @param df.expr.mut genes expression profiles with rows as observations and
#' columns as variables (nodes) of the tumor population
#' @param bootstrap if TRUE, uses a bootstrap or subsample routine to
#' compute the effects
#' @param runs bootstrap/subsampling runs
#' @param replace if TRUE, classic bootstrap, if FALSE subsampling without
#' replacement
#' @param frac the fraction of the data to sample, can either be one value
#' or a vector of two for different fractions for control and tumor,
#' respectively
#' populations
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
compute_differential_causal_effects <- function(graph.ctrl, df.expr.ctrl,
                                                graph.mut, df.expr.mut,
                                                bootstrap = FALSE, runs = 100,
                                                replace = FALSE, frac = 0.5) {
    if (bootstrap) {
        if (length(frac) == 1) { frac <- c(frac, frac) }
        ce.ctrl <- ce.mut <- 0
        for (b in seq_len(runs)) {
            df.expr.ctrl.sub <- df.expr.ctrl[
                sample(seq_len(nrow(df.expr.ctrl)),
                       ceiling(nrow(df.expr.ctrl)*frac[1]),
                       replace = replace), ]
            df.expr.mut.sub <- df.expr.mut[
                sample(seq_len(nrow(df.expr.mut)),
                       ceiling(nrow(df.expr.mut)*frac[2]),
                       replace = replace), ]
            ce.ctrl <- ce.ctrl+compute_causal_effects(graph.ctrl,
                                                      df.expr.ctrl.sub)
            ce.mut <- ce.mut+compute_causal_effects(graph.mut,
                                                    df.expr.mut.sub)
        }
        ce.ctrl <- ce.ctrl/runs
        ce.mut <- ce.mut/runs
    } else {
        ce.ctrl <- compute_causal_effects(graph.ctrl, df.expr.ctrl)
        ce.mut <- compute_causal_effects(graph.mut, df.expr.mut)
    }
    res <- t(as.matrix(ce.ctrl - ce.mut))
    class(res) <- "dce"
    return(res)
}

plot.dce <- function(x, dec=3, ...) {
    efreq <- efreqscale <- round(t(x)[which(t(x) != 0)], dec)
    efreqscale[abs(efreqscale) > 2] <- 2
    mnem::plotDnf(dnf, labels = efreq,
                  edgecol = rgb(abs(efreqscale)/2,0,(2-abs(efreqscale))/2))
}
