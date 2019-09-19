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
#' @author Kim Jablonski & Martin Pirkl
#' @return vector of causal effects
#' @export
#' @importFrom purrr map_dfc
#' @importFrom pcalg idaFast
#' @importFrom assertthat are_equal
#' @import graph tidyverse
#' @examples
compute_differential_causal_effects <- function(graph.ctrl, df.expr.ctrl,
                                                graph.mut, df.expr.mut,
                                                bootstrap = FALSE, runs = 100,
                                                replace = FALSE, frac = 0.5,
                                                strap = 0,
                                                method = "", ...) {
    if (bootstrap) {
        if (method %in% "full") {
            strap <- 1
        }
        if (length(frac) == 1) { frac <- c(frac, frac) }
        if (ceiling(nrow(df.expr.ctrl)*frac[1]) < ncol(df.expr.ctrl)) {
            frac[1] <- ncol(df.expr.ctrl)/nrow(df.expr.ctrl)
            print("too few control samples for subsampling at this fraction; reset to:")
            print(frac[1])
        }
        if (ceiling(nrow(df.expr.mut)*frac[2]) < ncol(df.expr.mut)) {
            frac[2] <- ncol(df.expr.mut)/nrow(df.expr.mut)
            print("too few tumor samples for subsampling at this fraction; reset to:")
            print(frac[2])
        }
        if (ncol(df.expr.mut) == nrow(df.expr.mut)) {
            print("too few samples for sampling with replacement")
            replace <- FALSE
        }
        if (strap) {
            dces <- 0
            for (b in seq_len(runs)) {
                if (!(method %in% "full")) {
                    df.expr.ctrl.sub <- df.expr.ctrl[
                        sample(seq_len(nrow(df.expr.ctrl)),
                               ceiling(nrow(df.expr.ctrl)*frac[1]),
                               replace = replace), ]
                    df.expr.mut.sub <- df.expr.mut[
                        sample(seq_len(nrow(df.expr.mut)),
                               ceiling(nrow(df.expr.mut)*frac[2]),
                               replace = replace), ]
                    dces <- dces +
                        compute_differential_causal_effects(graph.ctrl,
                                                            df.expr.ctrl.sub,
                                                            graph.mut,
                                                            df.expr.mut.sub)$dce
                } else {
                    dces <- dces + fulllin(graph.ctrl, df.expr.ctrl,
                                    graph.mut, df.expr.mut,
                                    ...)$dce
                }
            }
            res <- dces/runs
        } else {
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
            res <- t(as.matrix(ce.ctrl - ce.mut))
        }
    } else {
        if (!(method %in% "full")) {
            ce.ctrl <- compute_causal_effects(graph.ctrl, df.expr.ctrl)
            ce.mut <- compute_causal_effects(graph.mut, df.expr.mut)
            res <- t(as.matrix(ce.ctrl - ce.mut))
        } else {
            res <- fulllin(graph.ctrl, df.expr.ctrl,
                                    graph.mut, df.expr.mut,
                                    ...)$dce
        }
    }
    gtc <- as(graph.ctrl, "matrix")
    gtc[which(gtc != 0)] <- 1
    gtc <- mnem:::mytc(gtc)
    diag(gtc) <- 0
    res <- list(dce = res*gtc, graph = graph.ctrl, dcefull = res)
    class(res) <- "dce"
    return(res)
}
#' Plot dce object
#'
#' This function takes a differential causal effects object and plots
#' the dag with the dces
#' @param x dce object
#' @param dec rounding to dec decimals
#' @author Martin Pirkl
#' @method plot dce
#' @return plot of dag and dces
#' @export
#' @examples
plot.dce <- function(x, dec=3, ...) {
    x <- x$dce
    if(is.null(colnames(x))) {
        colnames(x) <- rownames(x) <- seq_len(ncol(x))
    }
    adj <- x
    adj[which(adj != 0)] <- 1
    adj <- mnem:::mytc(adj)
    diag(adj) <- 0
    dnf <- mnem:::adj2dnf(adj)
    dnf <- dnf[grep("=", dnf)]
    efreq <- efreqscale <- round(t(x)[which(t(x) != 0)], dec)
    efreqscale[abs(efreqscale) > 2] <- 2
    mnem::plotDnf(dnf, labels = efreq,
                  edgecol = rgb(abs(efreqscale)/2,0,(2-abs(efreqscale))/2))
}
fulllin <- function(g1, d1, g2, d2, conf = TRUE, ...) {
    mat1 <- as(g1, "matrix")
    mat2 <- as(g2, "matrix")
    mat1[which(mat1 != 0)] <- 1
    mat2[which(mat2 != 0)] <- 1
    dagtc <- mnem:::mytc(mat1)
    df <- rbind(d1, d2)
    colnames(df) <- paste0("X", seq_len(ncol(df)))
    df <- as.data.frame(cbind(df,
                              N = c(rep(1, nrow(d1)),
                                    rep(0, nrow(d2))),
                              T = c(rep(0, nrow(d1)),
                                    rep(1, nrow(d2)))))
    dce <- mat1*0
    for (i in seq_len(n)) {
        for (j in seq_len(n)) {
            if (dagtc[i, j] == 1 & i != j) {
                Z <- pcalg::backdoor(mat1, i, j, type = "dag")
                Z <- colnames(df)[Z]
                X <- colnames(df)[i]
                Y <- colnames(df)[j]
                if (length(Z) > 0 & conf) {
                    Lfit <- lm(paste0(Y, " ~ ",
                                      X, "*N + ",
                                      Z, "*N"
                                      ),
                               df)
                } else {
                    Lfit <- lm(paste0(Y, " ~ ",
                                      X, "*N"
                                      ),
                               df)
                }
                dce[i, j] <- Lfit$coefficients[grep(paste0(X, ":N"),
                                                    names(Lfit$coefficients))]
            }
        }
    }
    gtc <- as(g1, "matrix")
    gtc[which(gtc != 0)] <- 1
    gtc <- mnem:::mytc(gtc)
    diag(gtc) <- 0
    res <- list(dce = dce*gtc, graph = g1, dcefull = dce)
    class(res) <- "dce"
    return(res)
}


