test_that("something works", {
  set.seed(42)

  g <- as(matrix(c(0,0,1,0), 2, 2), "graphNEL")
  d <- simulate_data(g)

  g2 <- as(matrix(c(0,0,2.3,0), 2, 2), "graphNEL")
  d2 <- simulate_data(g2)

  res <- compute_differential_causal_effects(g, d, g2, d2)
  expect_equal(as.vector(res$dce), c(0, 0, -1.3, 0), tolerance=0.1)
})
