test_that("union works", {
    # create individual graphs
    g1 <- create_graph_from_dataframe(data.frame(
      from = c("A"),
      to = c("B")
    ), edge_weight = function() 0.5)
    g2 <- create_graph_from_dataframe(data.frame(
      from = c("A"),
      to = c("C")
    ), edge_weight = function() 1)
    g3 <- create_graph_from_dataframe(data.frame(
      from = c("X"),
      to = c("Y")
    ), edge_weight = function() -4)
    g4 <- create_graph_from_dataframe(data.frame(
      from = c("Y"),
      to = c("A")
    ), edge_weight = function() 6)

    # compute union
    g_union <- graph_union(c(g1, g2, g3, g4))

    # check result
    mat <- as(g_union, "matrix")

    expect_equal(dim(mat), c(5, 5))
    expect_equal(
        mat,
        matrix(
            c(
                0, 0, 0, 0, 6,
                0.5, 0, 0, 0, 0,
                1, 0, 0, 0, 0,
                0, 0, 0, 0, 0,
                0, 0, 0, -4, 0
            ),
            5, 5,
            dimnames = list(
                c("A", "B", "C", "X", "Y"),
                c("A", "B", "C", "X", "Y")
            )
        )
    )
})


test_that("DAG sampling works", {
  set.seed(42)

  # initial graph
  graph <- dce::create_random_DAG(100, prob = 0.5, eff_min = -4, eff_max = 3)
  mat <- as(graph, "matrix")

  expect_equal(dim(mat), c(100, 100))
  for (w in mat) {
    expect_lte(w, 3)
    expect_gte(w, -4)
  }

  # resampled graph
  graph_s <- dce::resample_edge_weights(graph, tp = 1, mineff = 4, maxeff = 7)
  mat_s <- as(graph_s, "matrix")

  expect_equal(dim(mat_s), c(100, 100))
  for (w in mat_s) {
    expect_lte(w, 7)
    expect_gte(w, 0)
  }
})
