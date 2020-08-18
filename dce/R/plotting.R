#' Plot network adjacency matrix
#'
#' Generic function which plots any adjacency matrix (assumes DAG)
#' @param adja.matrix Adjacency matrix of network
#' @param nodename.map node names
#' @param edge.colorscale.limits Limits for scale_edge_color_gradient2 (should contain 0). Useful to make plot comparable to others
#' @param nodesize Node sizes
#' @param labelsize Node label sizes
#' @param show.edge.labels Whether to show edge labels (DCEs)
#' @param use.symlog Scale edge colors using dce::symlog
#' @param highlighted.nodes List of nodes to highlight
#' @param legend.title Title of edge weight legend
#' @param value.matrix Optional matrix of edge weights if different from adjacency matrix
#' @param ... additional parameters
#' @author Martin Pirkl, Kim Philipp Jablonski
#' @return plot of dag and dces
#' @export
#' @import tidyverse ggraph purrr
#' @importFrom glue glue
#' @importFrom ggplot2 aes theme element_rect arrow unit coord_fixed scale_fill_manual
#' @importFrom tidygraph as_tbl_graph activate mutate
#' @importFrom rlang .data
plot_network <- function(
    adja.matrix,
    nodename.map = NULL, edgescale.limits = NULL,
    nodesize = 17, labelsize = 3,
    show.edge.labels = FALSE, use.symlog = FALSE,
    highlighted.nodes = c(),
    legend.title = "edge weight",
    value.matrix = NULL,
    ...
) {
    if (is.null(value.matrix)) {
        value.matrix <- adja.matrix
    }

    coords.dot <- purrr::map_dfr(
        Rgraphviz::agopen(as(adja.matrix, "graphNEL"), name="foo", layoutType="dot")@AgNode,
        function(node) {
            data.frame(x=node@center@x, y=node@center@y)
        }
    )

    if(is.null(edgescale.limits)) {
        edgescale.limits <- c(-max(abs(value.matrix), na.rm=TRUE), max(abs(value.matrix), na.rm=TRUE))
    }

    custom_breaks <- function(limits) {
        eps <- min(abs(limits)) / 100 # without this offset the outer breaks are somtetimes not shown
        breaks <- seq(limits[[1]] + eps, limits[[2]] - eps, length.out = 5)

        if(min(abs(limits)) < 1) {
            return(breaks)
        } else {
            return(round(breaks, 1))
        }
    }

    as_tbl_graph(as(adja.matrix, "graphNEL")) %>%
        activate(nodes) %>%
        mutate(
            label=if(is.null(nodename.map)) .data$name else nodename.map[.data$name],
            nodesize=nodesize,
            is.highlighted=.data$label %in% highlighted.nodes
        ) %>%
        activate(edges) %>%
        mutate(
            dce=pmap_dbl(
                list(.data$from, .data$to),
                function (f, t) { value.matrix[f, t] }
            ),
            dce.symlog=symlog(dce),
            label=.data$dce %>% round(2) %>% as.character
        ) %>%
    ggraph(layout=coords.dot) + # "sugiyama"
        geom_edge_diagonal(
            aes(
                color=if(use.symlog) .data$dce.symlog else .data$dce,
                alpha=abs(.data$dce),
                width=abs(.data$dce),
                label=if(show.edge.labels) .data$label else NULL,
                # start_cap = circle(.data$node1.nodesize, unit="native"), # uncomment once https://github.com/thomasp85/ggraph/pull/246 has been merged
                end_cap = circle(.data$node2.nodesize, unit="native")
            ),
            strength=0.5,
            arrow=arrow(type="closed", length=unit(3, "mm"))
        ) +
        geom_node_circle(aes(r=.data$nodesize, fill=.data$is.highlighted), color="black") +
        geom_node_text(aes(label=.data$label), size=labelsize) +
        coord_fixed() +
        scale_fill_manual(values=c("FALSE" = "white", "TRUE" = "red"), guide=FALSE) +
        scale_edge_color_gradient2(
            low="red", mid="grey", high="blue", midpoint=0,
            limits=if(use.symlog) symlog(edgescale.limits) else edgescale.limits,
            breaks=custom_breaks,
            name=if(use.symlog) glue("{legend.title} (symlog)") else legend.title
        ) +
        scale_edge_width(range=c(1, 3), limits=c(0, edgescale.limits[[2]]), guide=FALSE) +
        scale_edge_alpha(range=c(.1, 1), limits=c(0, edgescale.limits[[2]]), guide=FALSE) +
        theme(
            panel.background=element_rect(fill="white"),
            # legend.position="none"
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
    plot_network(x$graph, value.matrix = x$dce, legend.title = "DCE", ...)
}


#' @export
symlog <- function(x, base = 10, threshold = 1) {
    if(is.null(x)) {
        return(NULL)
    }

    ifelse(
        abs(x) < threshold,
        x,
        sign(x) * (threshold + log(sign(x) * x / threshold, base))
    )
}
