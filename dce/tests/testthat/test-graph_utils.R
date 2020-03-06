test_that("union works", {
    # create individual graphs
    g1 <- create_graph_from_dataframe(data.frame(
      from=c("A"),
      to=c("B")
    ), edge.weight=function() { 0.5 })
    g2 <- create_graph_from_dataframe(data.frame(
      from=c("A"),
      to=c("C")
    ), edge.weight=function() { 1 })
    g3 <- create_graph_from_dataframe(data.frame(
      from=c("X"),
      to=c("Y")
    ), edge.weight=function() { -4 })
    g4 <- create_graph_from_dataframe(data.frame(
      from=c("Y"),
      to=c("A")
    ), edge.weight=function() { 6 })

    # compute union
    g.union <- graph_union(c(g1, g2, g3, g4))

    # check result
    mat <- as(g.union, "matrix")

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
            dimnames=list(
                c("A", "B", "C", "X", "Y"),
                c("A", "B", "C", "X", "Y")
            )
        )
    )
})


test_that("DAG sampling works", {
  set.seed(42)

  # initial graph
  graph <- dce::create_random_DAG(100, prob = 0.5, lB = c(-7, -4), uB = c(2, 3))
  mat <- as(graph, "matrix")

  expect_equal(dim(mat), c(100, 100))
  for (w in mat[mat < 0]) {
    expect_lte(w, -4)
    expect_gte(w, -7)
  }
  for (w in mat[mat > 0]) {
    expect_lte(w, 3)
    expect_gte(w, 2)
  }

  # resampled graph
  graph.s <- dce::resample_edge_weights(graph, lB = c(-10, -9), uB = c(6, 8))
  mat.s <- as(graph.s, "matrix")

  expect_equal(dim(mat.s), c(100, 100))
  for (w in mat.s[mat.s < 0]) {
    expect_lte(w, -9)
    expect_gte(w, -10)
  }
  for (w in mat.s[mat.s > 0]) {
    expect_lte(w, 8)
    expect_gte(w, 6)
  }
})
