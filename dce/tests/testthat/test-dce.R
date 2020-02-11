test_that("positive beta can be recovered", {
  set.seed(42)

  graph.wt <- as(matrix(c(0, 0, 1e-42, 0), 2, 2), "graphNEL")
  X.wt <- simulate_data(graph.wt)

  graph.mt <- as(matrix(c(0, 0, 1.5, 0), 2, 2), "graphNEL")
  X.mt <- simulate_data(graph.mt)

  res <- compute_differential_causal_effects(graph.wt, X.wt, graph.mt, X.mt)

  expect_equal(as.vector(res$dce), c(0, 0, 1.5, 0), tolerance=0.1)
})


test_that("negative beta can be recovered", {
  set.seed(42)

  graph.wt <- as(matrix(c(0, 0, 1e-42, 0), 2, 2), "graphNEL")
  X.wt <- simulate_data(graph.wt)

  graph.mt <- igraph::igraph.to.graphNEL(igraph::graph_from_adjacency_matrix(matrix(c(0, 0, -1.5, 0), 2, 2)))
  X.mt <- simulate_data(graph.mt)

  res <- compute_differential_causal_effects(graph.wt, X.wt, graph.mt, X.mt)

  expect_equal(as.vector(res$dce), c(0, 0, -1.5, 0), tolerance=0.1)
})
