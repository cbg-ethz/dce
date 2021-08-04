test_that("different p-value computations work", {
  set.seed(42)

  node_names <- c("A", "B", "C")

  graph_wt <- matrix(c(0, 0, 0, 1e-42, 0, 0, 1e-42, 0, 0), 3, 3)
  rownames(graph_wt) <- colnames(graph_wt) <- node_names
  X_wt <- simulate_data(graph_wt)

  graph_mt <- matrix(c(0, 0, 0, 1.5, 0, 0, -.7, 0, 0), 3, 3)
  rownames(graph_mt) <- colnames(graph_mt) <- node_names
  X_mt <- simulate_data(graph_mt)

  res_lr <- dce::dce(graph_wt, X_wt, X_mt, test = "lr")
  expect_lt(res_lr$dce_pvalue["A", "B"], 1e-100)

  res_wald <- dce::dce(graph_wt, X_wt, X_mt, test = "wald")
  expect_lt(res_wald$dce_pvalue["A", "B"], 1e-100)

  res_vcov <- dce::dce(graph_wt, X_wt, X_mt, test = "vcovHC")
  expect_lt(res_vcov$dce_pvalue["A", "B"], 1e-100)
})
