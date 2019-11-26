#' Compute the true casual effects of a simulated dag
#'
#' This function takes a DAG with edgeweights as input and computes
#' the causal effects of all nodes on all direct and indirect children in the
#' DAG. Alternatively see pcalg::causalEffect for pairwise computation.
#' @param g graphNEL object
#' @author Martin Pirkl
#' @return matrix of causal effects
#' @export
#' @importFrom pcalg causalEffect
#' @import graph tidyverse
#' @importFrom expm %^%
#' @examples
#' graph.wt <- as(matrix(c(0,0,0,1,0,0,0,1,0), 3), "graphNEL")
#' trueEffects(graph.wt)
trueEffects <- function(g) {
    ## n <- ncol(as(g, "matrix"))
    ## te <- matrix(0, n, n)
    ## for (i in seq_len(n-1)) {
    ##     for (j in i:n) {
    ##         te[i, j] <- causalEffect(g, j, i)
    ##     }
    ## }
    ## return(te)
    a <- as(g, "matrix")
    ae <- a
    for (i in 2:nrow(a)) {
        ae <- ae + a%^%i
    }
    return(ae)
}

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
#' @import graph tidyverse stats
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
    ce.t <- purrr::map_dfc(seq_len(n), function(x) {
        pcalg::idaFast(x, seq_len(n), cov.mat, graph)
    })

    t(ce.t)
}
#' Compute the differential causal effects
#'
#' This function takes a DAG and gene expression data of two distinct types
#' of populations (e.g. tumor versus normal) as input and computes
#' the differential causal effects of all nodes on all other nodes in the DAG
#' between the tow populations of samples.
#' @param graph.ctrl DAG as a graphNEL object of the control population
#' @param df.expr.ctrl Gene expression profiles with rows as observations and
#' columns as variables (nodes) of the control population
#' @param graph.mut DAG as a graphNEL object of the tumor population
#' @param df.expr.mut genes expression profiles with rows as observations and
#' columns as variables (nodes) of the tumor population
#' @param method Either "full" for the full linear model using all
#'     samples together to directly compute the differential causal effect
#'     or "part" for two linear models on both conditions to first compute
#'     the causal effect and then the difference
#' @param bootstrap if TRUE, uses a bootstrap or subsample routine to
#' compute the effects
#' @param runs bootstrap/subsampling runs
#' @param replace if TRUE, classic bootstrap, if FALSE subsampling without
#' replacement
#' @param frac the fraction of the data to sample, can either be one value
#' or a vector of two for different fractions for control and tumor,
#' respectively populations
#' @param bootMethod if "diff" (default) bootstraps the differential causal
#' effects, if "cause", bootstraps on the causal effects
#' (not possible for "full")
#' @param ... further arguments passed to `fulllin`
#' @author Kim Jablonski & Martin Pirkl
#' @return vector of causal effects
#' @export
#' @importFrom purrr map_dfc
#' @importFrom pcalg idaFast
#' @importFrom assertthat are_equal
#' @importFrom methods as
#' @import graph tidyverse
#' @examples
#' graph.wt <- as(matrix(c(0,0,0,1,0,0,0,1,0), 3), "graphNEL")
#' graph.mt <- resample_edge_weights(graph.wt)
#' X.wt <- simulate_data(graph.wt)
#' X.mt <- simulate_data(graph.mt)
#' compute_differential_causal_effects(graph.wt, X.wt, graph.mt, X.mt)
compute_differential_causal_effects <- function(
    graph.ctrl, df.expr.ctrl,
    graph.mut, df.expr.mut,
    method = "full",
    bootstrap = FALSE, runs = 100,
    replace = FALSE, frac = 0.5,
    bootMethod = "diff", ...
) {
    if (bootstrap) {
        if (method %in% "full") {
            bootMethod <- "diff"
        }
        if (length(frac) == 1) { frac <- c(frac, frac) }
        if (ceiling(nrow(df.expr.ctrl)*frac[1]) < ncol(df.expr.ctrl)) {
            frac[1] <- ncol(df.expr.ctrl)/nrow(df.expr.ctrl)
            print(
                paste0("too few control samples for subsampling ",
                "at this fraction; reset to:")
            )
            print(frac[1])
        }
        if (ceiling(nrow(df.expr.mut)*frac[2]) < ncol(df.expr.mut)) {
            frac[2] <- ncol(df.expr.mut)/nrow(df.expr.mut)
            print(
                paste0("too few tumor samples for subsampling ",
                "at this fraction; reset to:")
            )
            print(frac[2])
        }
        if (ncol(df.expr.mut) == nrow(df.expr.mut)) {
            print("too few samples for sampling with replacement")
            replace <- FALSE
        }
        if (bootMethod %in% "diff") {
            dces <- 0
            for (b in seq_len(runs)) {
                df.expr.ctrl.sub <- df.expr.ctrl[
                    sample(
                        seq_len(nrow(df.expr.ctrl)),
                        ceiling(nrow(df.expr.ctrl)*frac[1]),
                        replace = replace
                    ),
                ]
                df.expr.mut.sub <- df.expr.mut[
                    sample(
                        seq_len(nrow(df.expr.mut)),
                        ceiling(nrow(df.expr.mut)*frac[2]),
                        replace = replace
                    ),
                ]
                if (!(method %in% "full")) {
                    dces <- dces +
                        compute_differential_causal_effects(
                            graph.ctrl,
                            df.expr.ctrl.sub,
                            graph.mut,
                            df.expr.mut.sub
                        )$dce
                } else {
                    dces <- dces + fulllin(
                        graph.ctrl, df.expr.ctrl.sub,
                        graph.mut, df.expr.mut.sub,
                        ...
                    )$dce
                }
            }
            res <- dces/runs
        } else {
            ce.ctrl <- ce.mut <- 0
            for (b in seq_len(runs)) {
                df.expr.ctrl.sub <- df.expr.ctrl[
                    sample(
                        seq_len(nrow(df.expr.ctrl)),
                        ceiling(nrow(df.expr.ctrl)*frac[1]),
                        replace = replace
                    ),
                ]
                df.expr.mut.sub <- df.expr.mut[
                    sample(
                        seq_len(nrow(df.expr.mut)),
                        ceiling(nrow(df.expr.mut)*frac[2]),
                        replace = replace
                    ),
                ]
                ce.ctrl <- ce.ctrl+compute_causal_effects(
                    graph.ctrl,
                    df.expr.ctrl.sub
                )
                ce.mut <- ce.mut+compute_causal_effects(
                    graph.mut,
                    df.expr.mut.sub
                )
            }
            ce.ctrl <- ce.ctrl/runs
            ce.mut <- ce.mut/runs
            res <- as.matrix(ce.ctrl - ce.mut)
        }
    } else {
        if (!(method %in% "full")) {
            ce.ctrl <- compute_causal_effects(graph.ctrl, df.expr.ctrl)
            ce.mut <- compute_causal_effects(graph.mut, df.expr.mut)
            res <- as.matrix(ce.ctrl - ce.mut)
        } else {
            res <- fulllin(
                graph.ctrl, df.expr.ctrl,
                graph.mut, df.expr.mut,
                ...
            )$dce
        }
    }
    res <- list(dce = res, graph = graph.ctrl)

    rownames(res$dce) <- nodes(graph.ctrl)
    colnames(res$dce) <- nodes(graph.ctrl)

    class(res) <- "dce"
    return(res)
}
#' Plot dce object
#'
#' This function takes a differential causal effects object and plots
#' the dag with the dces
#' @param x dce object
#' @param ... additional parameters
#' @author Martin Pirkl, Kim Philipp Jablonski
#' @method plot dce
#' @return plot of dag and dces
#' @export
#' @import tidyverse ggraph purrr
#' @importFrom ggplot2 aes theme element_rect arrow unit
#' @importFrom tidygraph as_tbl_graph activate mutate
#' @importFrom rlang .data
plot.dce <- function(x, ...) {
    as_tbl_graph(x$graph) %>%
        activate(edges) %>%
        mutate(
            dce=pmap_dbl(
                list(.data$from, .data$to),
                function (f, t) { x$dce[f, t] }
            ),
            label=.data$dce %>% round(2) %>% as.character
        ) %>%
    ggraph(layout="sugiyama") +
        geom_edge_diagonal(
            aes(
                label=.data$label, width=abs(.data$dce), color=.data$dce,
                start_cap=label_rect(.data$node1.name),
                end_cap=label_rect(.data$node2.name)
            ),
            strength=0.5,
            arrow=arrow(length=unit(3, "mm"))
        ) +
        geom_node_point(color="grey", size=8) +
        geom_node_text(aes(label=.data$name)) +
        scale_edge_color_gradient2(
            low="red", mid="grey", high="blue",
            midpoint=0
        ) +
        scale_edge_width(range=c(1, 3)) +
        theme(
            panel.background=element_rect(fill="white"),
            legend.position="none"
        )
}
