#' Compute statistic for permutation test
#'
#' This function takes the same parameters as
#' compute_differential_effects in addition to
#' control parameters for the permutations.
#' @param graph.ctrl DAG as a graphNEL object of the control population
#' @param df.expr.ctrl Gene expression profiles with rows as observations and
#' columns as variables (nodes) of the control population
#' @param graph.mut DAG as a graphNEL object of the tumor population
#' @param df.expr.mut genes expression profiles with rows as observations and
#' columns as variables (nodes) of the tumor populationt
#' @param runs number of permutations
#' @param statistic function, which computes a single valued
#' statistic for each run
#' @param ... additional parameters for compute_differential_effects
#' @author Martin Pirkl
#' @return numeric vector with permutation statistics of length runs
#' @export
#' @examples
#' 1
compute_permutations <- function(normal, dn, tumor, dt, runs=10,
                                 statistic = function(x)
                                     return(sum(abs(x))),
                                 ...) {
    statistics <- numeric(runs)
    stats_edges <- list()
    for (i in seq_len(runs)) {
        dnp <- dn
        colnames(dnp) <- sample(colnames(dn), ncol(dn))
        dtp <- dt
        colnames(dtp) <- sample(colnames(dt), ncol(dt))
        dnp <- dnp[, order(colnames(dnp))]
        dtp <- dtp[, order(colnames(dtp))]
        dceip <- compute_differential_causal_effects(normal, dnp,
                                                     tumor, dtp, ...
                                                     )
        statistics[i] <- statistic(dceip$dce)
        stats_edges[[i]] <- as.vector(dceip$dce)
    }
    return(list(statistics = statistics, stats_edges = stats_edges))
}
#' Compute the true casual effects of a simulated dag
#'
#' This function takes a DAG with edgeweights as input and computes
#' the causal effects of all nodes on all direct and indirect children in the
#' DAG. Alternatively see pcalg::causalEffect for pairwise computation.
#' @param g graphNEL object
#' @param partial if FALSE computes the total causal effects and not just
#' the partial edge effects
#' @author Martin Pirkl
#' @return matrix of causal effects
#' @export
#' @importFrom pcalg causalEffect
#' @import graph tidyverse
#' @importFrom expm %^%
#' @examples
#' graph.wt <- as(matrix(c(0,0,0,1,0,0,0,1,0), 3), "graphNEL")
#' trueEffects(graph.wt)
trueEffects <- function(g, partial = FALSE) {
    ## n <- ncol(as(g, "matrix"))
    ## te <- matrix(0, n, n)
    ## for (i in seq_len(n-1)) {
    ##     for (j in i:n) {
    ##         te[i, j] <- causalEffect(g, j, i)
    ##     }
    ## }
    ## return(te)
    a <- as(g, "matrix")
    if (partial) {
        ae <- a
    } else {
        ae <- a
        for (i in 2:nrow(a)) {
            ae <- ae + a%^%i
        }
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
#' @param theta prior estimated theta value
#' @param partial if TRUE only computes the partial causal effects on
#' the edges, else computes the total causal effect
#' @param conf if TRUE accounts for confounders
#' @param link.log.base base of logarithm of link
#' @param ... additional parameters for glm2::glm2
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
    bootMethod = "diff", errDist = "nbinom",
    theta = NULL, partial = FALSE,
    link.log.base = exp(1),
    ...
) {
    if (is.null(theta)) {
        theta <- estimateTheta(rbind(df.expr.ctrl, df.expr.mut))
    }
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
                        errDist = errDist, theta = theta,
                        partial = partial,
                        link.log.base = link.log.base,
                        ...
                    )$dce
                }
            }
            res <- dces/runs
            res.p <- NULL
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
            res.p <- NULL
        }
    } else {
        if (!(method %in% "full")) {
            ce.ctrl <- compute_causal_effects(graph.ctrl, df.expr.ctrl)
            ce.mut <- compute_causal_effects(graph.mut, df.expr.mut)
            res <- as.matrix(ce.ctrl - ce.mut)
            res.p <- NULL
        } else {
            tmp <- fulllin(
                graph.ctrl, df.expr.ctrl,
                graph.mut, df.expr.mut,
                errDist = errDist, theta = theta,
                partial = partial,
                link.log.base = link.log.base,
                ...
            )

            res <- tmp$dce
            res.p <- tmp$dce.p
        }
    }
    mat <- as(graph.ctrl, "matrix")
    mat[which(mat != 0)] <- 1
    dagtc <- nem::transitive.closure(mat, mat=TRUE)
    res <- res*dagtc
    res.p <- res.p*dagtc
    out <- list(dce = res, dce.p = res.p, graph = graph.ctrl,
                theta = theta)
    rownames(out$dce) <- nodes(graph.ctrl)
    colnames(out$dce) <- nodes(graph.ctrl)
    class(out) <- "dce"
    return(out)
}
#' Plot dce object
#'
#' This function takes a differential causal effects object and plots
#' the dag with the dces
#' @param x dce object
#' @param nodename.map node names
#' @param ... additional parameters
#' @author Martin Pirkl, Kim Philipp Jablonski
#' @method plot dce
#' @return plot of dag and dces
#' @export
#' @import tidyverse ggraph purrr
#' @importFrom ggplot2 aes theme element_rect arrow unit
#' @importFrom tidygraph as_tbl_graph activate mutate
#' @importFrom rlang .data
plot.dce <- function(x, nodename.map = NULL, ...) {
    as_tbl_graph(x$graph) %>%
        activate(nodes) %>%
        mutate(
            label=if(is.null(nodename.map)) .data$name else nodename.map[.data$name]
        ) %>%
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
                alpha=abs(.data$dce),
                start_cap=label_rect(.data$node1.name),
                end_cap=label_rect(.data$node2.name)
            ),
            strength=0.5,
            arrow=arrow(length=unit(3, "mm"))
        ) +
        geom_node_point(color="grey", size=8) +
        geom_node_text(aes(label=.data$label)) +
        scale_edge_color_gradient2(
            low="red", mid="violet", high="blue",
            midpoint=0
        ) +
        scale_edge_width(range=c(1, 3)) +
        scale_edge_alpha(range=c(.1, 1)) +
        theme(
            panel.background=element_rect(fill="white"),
            legend.position="none"
        )
}
#' Compute pathway enrichment
#'
#' This function computes a p-value.
#' @param graph DAG as a graphNEL object
#' @param X.wt Expression values of wild type (if dce=NULL)
#' @param X.mt Expression values of mutant (if dce=NULL)
#' @param statistic Statistic to compute
#' @param permutation_count How many permutations to do
#' @param pvalue.method "hmp" for harmonic median, "perm" for permutation test
#' @param theta prior estimated theta value
#' @param partial if TRUE only computes the partial causal
#' effects (i.e. edge weights)
#' @param dce optional dce object
#' the edges, else computes the total causal effect
#' @param ... additional parameters for compute_differential_causal_effects
#' @author Hinz und Kunz
#' @return Enrichment p-value
#' @export
#' @examples
#' graph.wt <- as(matrix(c(0,0,0,1,0,0,0,1,0), 3), "graphNEL")
#' graph.mt <- resample_edge_weights(graph.wt)
#' X.wt <- simulate_data(graph.wt)
#' X.mt <- simulate_data(graph.mt)
#' compute_enrichment(graph.wt, X.wt, X.mt)
compute_enrichment <- function(
    graph, X.wt, X.mt,
    statistic = function(x) { sum(abs(x)) },
    permutation_count = 100,
    pvalue.method = "hmp", # "perm"
    theta = NULL, partial = FALSE,
    dce = NULL, ...
) {
    if (is.null(dce)) {
        # compute observed statistic
        res <- compute_differential_causal_effects(
            graph, X.wt, graph, X.mt,
            theta = theta,
            partial = partial, ...
        )
    } else {
        res <- dce
    }

    if (pvalue.method == "hmp") {
        # aggregate p-values
        tmp <- res$dce.p[!is.na(res$dce.p)]
        p.val <- as.numeric(harmonicmeanp::p.hmp(tmp, L = length(tmp)))

        return(list(p.value = p.val, p.edges = res$dce.p))
    } else if (pvalue.method == "perm") {
        # compute empirical p-value
        dce.inferred <- res$dce
        stats.inferred <- statistic(dce.inferred)

        stats.permuted <- compute_permutations(
            graph, X.wt, graph, X.mt,
            runs = permutation_count,
            statistic = statistic,
            theta = theta,
            partial = partial,
            ...
        )

        p.value <- sum(stats.permuted[[1]] >= stats.inferred) / permutation_count

        ## for edges
        g.vec <- as.vector(as(graph, "matrix"))
        stats_mat <- rbind(as.vector(dce.inferred),
                           do.call("rbind", stats.permuted[[2]]))
        stats_mat <- stats_mat[, which(g.vec != 0)]
        p.edges <- apply(stats_mat, 2, function(x) {
            x <- abs(x)
            y <- sum(x[-1] >= x[1])/(length(x)-1)
            return(y)
        })
        g.vec[which(g.vec != 0)] <- p.edges
        g.vec[which(g.vec == 0)] <- NA
        return(list(p.value = p.value, p.edges = g.vec))
    }
}
