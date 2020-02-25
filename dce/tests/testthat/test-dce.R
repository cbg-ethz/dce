test_that("positive beta can be recovered", {
  set.seed(42)

  node.names <- c("A", "B")

  graph.wt <- matrix(c(0, 0, 1e-42, 0), 2, 2)
  rownames(graph.wt) <- colnames(graph.wt) <- node.names
  X.wt <- simulate_data(graph.wt)

  graph.mt <- matrix(c(0, 0, 1.5, 0), 2, 2)
  rownames(graph.mt) <- colnames(graph.mt) <- node.names
  X.mt <- simulate_data(graph.mt)

  res <- dce::dce.nb(graph.wt, X.wt, X.mt)
  res

  expect_equal(as.vector(res$dce), c(NA, NA, 1.5, NA), tolerance = 0.1)

  expect_equal(rownames(res$dce), node.names)
  expect_equal(colnames(res$dce), node.names)
  expect_equal(rownames(res$dce.pvalue), node.names)
  expect_equal(colnames(res$dce.pvalue), node.names)
})


test_that("negative beta can be recovered", {
  set.seed(42)

  graph.wt <- matrix(c(0, 0, 1e-42, 0), 2, 2)
  X.wt <- simulate_data(graph.wt)

  graph.mt <- matrix(c(0, 0, -1.5, 0), 2, 2)
  X.mt <- simulate_data(graph.mt)

  res <- dce::dce.nb(graph.wt, X.wt, X.mt)
  res

  expect_equal(as.vector(res$dce), c(NA, NA, -1.5, NA), tolerance = 0.1)
})


test_that("igraph input works", {
  set.seed(42)
  devtools::load_all("/Users/kimja/university/PhD/projects/causal_pathways/dce/")

  graph <- igraph::make_tree(3, children = 3, mode = "out")
  graph.wt <- igraph::set.edge.attribute(graph, name = "weight", value = 1)
  graph.mt <- igraph::set.edge.attribute(graph, name = "weight", value = 2.4)

  X.wt <- simulate_data(graph.wt, n = 1000)
  X.mt <- simulate_data(graph.mt, n = 1000)

  res <- dce::dce.nb(graph, X.wt, X.mt)
  res

  expect_equal(as.vector(res$dce), c(NA, NA, NA, 1.4, NA, NA, 1.4, NA, NA), tolerance = 0.1)
})


test_that("graphNEL input works", {
  set.seed(42)

  graph.wt <- as(matrix(c(0, 0, 1e-42, 0), 2, 2), "graphNEL")
  X.wt <- simulate_data(graph.wt)

  graph.mt <- as(matrix(c(0, 0, 1.6, 0), 2, 2), "graphNEL")
  X.mt <- simulate_data(graph.mt)

  res <- dce::dce.nb(graph.wt, X.wt, X.mt)
  res

  expect_equal(as.vector(res$dce), c(NA, NA, 1.6, NA), tolerance = 0.1)
})


test_that("graph nodes and simulated data column mismatch throws error", {
  set.seed(42)

  graph <- dce::create_graph_from_dataframe(data.frame(
    from=c("A"),
    to=c("B")
  ))

  X.wt <- simulate_data(graph)
  X.mt <- simulate_data(graph)

  colnames(X.wt) <- c("Hinz", "Kunz") # introduce error

  expect_error(
    dce::dce.nb(graph, X.wt, X.mt),
    "Not all nodes have expression vector in WT data"
  )
})


test_that("adjustment sets work", {
  set.seed(42)

  graph <- dce::create_graph_from_dataframe(data.frame(
    from=c("A", "B"),
    to=c("C", "C")
  ))

  expect_equal(
    get.adjustment.set(as(graph, "matrix"), which(nodes(graph) == "A"), which(nodes(graph) == "B"), "parents"),
    vector(mode="character")
  )
  expect_equal(
    get.adjustment.set(as(graph, "matrix"), which(nodes(graph) == "A"), which(nodes(graph) == "B"), "minimal"),
    c("C")
  )
})
