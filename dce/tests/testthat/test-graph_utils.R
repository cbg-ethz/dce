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
