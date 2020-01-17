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


#' @export
graph_union <- function(graph1, graph2) {
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
