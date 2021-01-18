test_that("simple propagation works", {
  set.seed(42)

  graph <- create_graph_from_dataframe(data.frame(
    from = c("A", "B"),
    to = c("B", "C")
  ), edge_weight = function() { 10 })

  X <- simulate_data(graph, n = 1000, dist_mean = 42)

  X_summary <- X %>%
    as.data.frame %>%
    summarise_all(list(mean))

  expect_equal(X_summary$A, 42, tolerance = 1e-2)
  expect_gt(X_summary$B, 42)
  expect_gt(X_summary$C, 42)
})
