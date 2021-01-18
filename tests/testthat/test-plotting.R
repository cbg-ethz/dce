test_that("trivial network plot works", {
  dce::plot_network(
    adja_matrix = matrix(c(0, 0, 0, 1, 0, 0, 0, 1, 0), 3, 3)
  )

  dce::plot_network(
    adja_matrix = matrix(c(0, 0, 0, 1, 0, 0, 0, 1, 0), 3, 3),
    value_matrix = matrix(c(0, 0, 0, 1.7, 0, 0, 0, NA, 0), 3, 3)
  )
})
