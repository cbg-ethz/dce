#' Plot network adjacency matrix
#'
#' Generic function which plots any adjacency matrix (assumes DAG)
#' @param adja_matrix Adjacency matrix of network
#' @param nodename_map node names
#' @param edge.colorscale.limits Limits for scale_edge_color_gradient2 (should contain 0). Useful to make plot comparable to others
#' @param nodesize Node sizes
#' @param labelsize Node label sizes
#' @param show_edge_labels Whether to show edge labels (DCEs)
#' @param use_symlog Scale edge colors using dce::symlog
#' @param highlighted_nodes List of nodes to highlight
#' @param legend_title Title of edge weight legend
#' @param value_matrix Optional matrix of edge weights if different from adjacency matrix
#' @param ... additional parameters
#' @author Martin Pirkl, Kim Philipp Jablonski
#' @return plot of dag and dces
#' @export
#' @import tidyverse ggraph purrr
#' @importFrom glue glue
#' @importFrom ggplot2 aes theme element_rect arrow unit coord_fixed scale_fill_manual waiver
#' @importFrom tidygraph as_tbl_graph activate mutate
#' @importFrom rlang .data
#' @importFrom igraph graph_from_adjacency_matrix
plot_network <- function(
    adja_matrix,
    nodename_map = NULL, edgescale_limits = NULL,
    nodesize = 17, labelsize = 3,
    show_edge_labels = FALSE, use_symlog = FALSE,
    highlighted_nodes = c(),
    legend_title = "edge weight",
    value_matrix = NULL,
    ...
) {
    # sanitize input
    if (is.null(value_matrix)) {
        value_matrix <- adja_matrix
    }

    if (is.null(rownames(adja_matrix)) || is.null(colnames(adja_matrix))) {
        warning(
            "Not nodenames set, using dummy names...",
            call. = FALSE
        )

        node_names <- seq_len(dim(adja_matrix)[[1]])
        rownames(adja_matrix) <- colnames(adja_matrix) <- node_names
    }

    # compute node coordinates
    tmp <- adja_matrix
    tmp[tmp != 0] <- 1

    coords_dot <- purrr::map_dfr(
        Rgraphviz::agopen(
            as(tmp, "graphNEL"),
            name = "foo",
            layoutType = "dot"
        )@AgNode,
        function(node) {
            data.frame(x = node@center@x, y = node@center@y)
        }
    )

    # handle scale setup
    if (is.null(edgescale_limits)) {
        edgescale_limits <- c(
            -max(abs(value_matrix), na.rm = TRUE),
            max(abs(value_matrix), na.rm = TRUE)
        )
    }

    custom_breaks <- function(limits) {
        if (any(is.infinite(limits))) {
            return(ggplot2::waiver())
        }

        eps <- min(abs(limits)) / 100 # without this offset the outer breaks are somtetimes not shown
        breaks <- seq(limits[[1]] + eps, limits[[2]] - eps, length.out = 5)

        if (min(abs(limits)) < 1) {
            return(breaks)
        } else {
            return(round(breaks, 1))
        }
    }

    # create plot
    as_tbl_graph(igraph::graph_from_adjacency_matrix(
        adja_matrix, weighted = TRUE
    )) %>%
        activate(nodes) %>%
        mutate(
            label = if (is.null(nodename_map)) .data$name else nodename_map[.data$name],
            nodesize = nodesize,
            is.highlighted = .data$label %in% highlighted_nodes
        ) %T>%
        with({
            label_list <- as.data.frame(.)$label
            extra_nodes <- setdiff(highlighted_nodes, label_list)

            if (length(extra_nodes) > 0) {
                label_str <- glue::glue_collapse(extra_nodes, sep = ", ")
                warning(
                    glue::glue("Invalid highlighted nodes: {label_str}"),
                    call. = FALSE
                )
            }
        }) %>%
        activate(edges) %>%
        mutate(
            dce = pmap_dbl(
                list(.data$from, .data$to),
                function(f, t) { value_matrix[f, t] }
            ),
            dce.symlog = symlog(dce),
            label = .data$dce %>% round(2) %>% as.character
        ) %>%
    ggraph(layout = coords_dot) + # "sugiyama"
        geom_node_circle(
            aes(r = .data$nodesize, fill = .data$is.highlighted),
            color = "black"
        ) +
        geom_edge_diagonal(
            aes(
                color = if (use_symlog) .data$dce.symlog else .data$dce,
                alpha = abs(.data$dce),
                width = abs(.data$dce),
                label = if (show_edge_labels) .data$label else NULL,
                linetype = is.na(.data$dce),
                # proper caps can be used again when https://github.com/thomasp85/ggraph/issues/254 is fixed
                # start_cap = circle(.data$node1.nodesize, unit="native"),
                # end_cap = circle(.data$node2.nodesize, unit="native")
            ),
            strength = 0.5,
            arrow = arrow(type = "closed", length = unit(3, "mm"))
        ) +
        geom_node_text(aes(label = .data$label), size = labelsize) +
        coord_fixed() +
        scale_fill_manual(
            values = c("FALSE" = "white", "TRUE" = "red"), guide = FALSE
        ) +
        scale_edge_color_gradient2(
            low = "red", mid = "grey", high = "blue", midpoint = 0,
            limits = if (use_symlog) symlog(edgescale_limits) else edgescale_limits,
            breaks = custom_breaks,
            name = if (use_symlog) glue("{legend_title} (symlog)") else legend_title,
            na.value = "black"
        ) +
        scale_edge_width(
            range = c(1, 3), limits = c(0, edgescale_limits[[2]]),
            na.value = 1,
            guide = FALSE
        ) +
        scale_edge_alpha(
            range = c(.1, 1), limits = c(0, edgescale_limits[[2]]),
            na.value = 1,
            guide = FALSE
        ) +
        scale_edge_linetype_manual(
            values = c("FALSE" = "solid", "TRUE" = "dashed"),
            guide = FALSE
        ) +
        theme(
            panel.background = element_rect(fill = "white"),
            # legend.position = "none"
        )
}

#' Plot dce object
#'
#' This function takes a differential causal effects object and plots
#' the dag with the dces
#' @param x dce object
#' @param ... Parameters passed to dce::plot_network
#' @author Martin Pirkl, Kim Philipp Jablonski
#' @method plot dce
#' @return plot of dag and dces
#' @export
plot.dce <- function(x, ...) {
    plot_network(x$graph, value_matrix = x$dce, legend_title = "DCE", ...)
}


#' @export
symlog <- function(x, base = 10, threshold = 1) {
    if (is.null(x)) {
        return(NULL)
    }

    ifelse(
        abs(x) < threshold,
        x,
        sign(x) * (threshold + log(sign(x) * x / threshold, base))
    )
}
