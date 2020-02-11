#' @export
create_graph_from_dataframe <- function(
  df.graph,
  edge.weight=function() { runif(1, 0.5, 2) }
) {
  g <- tidygraph::as_tbl_graph(df.graph)
  ig <- igraph::as.igraph(g)
  graph <- igraph::igraph.to.graphNEL(ig)

  w <- graph@edgeData@data
  for (i in 1:length(w)) {
    w[[i]]$weight <- edge.weight()
  }
  graph@edgeData@data <- w

  return(graph)
}


#' @noRd
graph_union_two <- function(graph1, graph2) {
  # create union
  graph.m <- igraph::union(
    igraph::igraph.from.graphNEL(graph1),
    igraph::igraph.from.graphNEL(graph2)
  )

  # merge edge weight attributes (if needed)
  if (is.null(igraph::get.edge.attribute(graph.m, "weight"))) {
    attr1 <- igraph::get.edge.attribute(graph.m, "weight_1")
    attr2 <- igraph::get.edge.attribute(graph.m, "weight_2")
    tmp <- igraph::set.edge.attribute(
      graph.m, "weight",
      value = ifelse(is.na(attr1), attr2, attr1)
    )

    # remove superfluous weight attributes (may cause errors on later conversions)
    tmp <- igraph::remove.edge.attribute(tmp, "weight_1")
    tmp <- igraph::remove.edge.attribute(tmp, "weight_2")
  } else {
    tmp <- graph.m
  }

  # return result
  gn <- igraph::igraph.to.graphNEL(tmp)
  return(gn)
}


#' @export
graph_union <- function(graph_list) {
  Reduce(graph_union_two, graph_list)
}

#' @export
propagate_gene_edges <- function(graph) {
  # propagate edges
  ig <- igraph::igraph.from.graphNEL(graph)
  vertex.names <- igraph::vertex_attr(ig, "name")

  for (source.idx in igraph::V(ig)) {
    source <- vertex.names[source.idx]

    for (target.idx in igraph::neighbors(ig, source.idx, mode="out")) {
      target <- vertex.names[target.idx]

      if (substr(target, start=0, stop=3) != "hsa") {
        # source is not connected to gene

        for (bridge.idx in igraph::neighbors(ig, target.idx, mode="out")) {
          bridge <- vertex.names[bridge.idx]
          stopifnot(substr(bridge, start=0, stop=3) == "hsa")

          ig <- igraph::add.edges(ig, c(source.idx, bridge.idx))
        }
      }
    }
  }

  graph.prop <- igraph::igraph.to.graphNEL(ig)

  # remove non-gene nodes
  hsa.nodes <- Filter(function(x) { substr(x, 0, 3) == "hsa" }, nodes(graph.prop))
  graph.filter <- graph::subGraph(hsa.nodes, graph.prop)

  return(graph.filter)
}
