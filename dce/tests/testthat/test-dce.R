test_that("something works", {
  g <- randomDAG(10, .5)
  d <- rmvDAG_2(100, g)

  g2 <- newWeights(g)
  d2 <- rmvDAG_2(100, g2)

  compute_differential_causal_effects(g, d, g2, d2)
})
