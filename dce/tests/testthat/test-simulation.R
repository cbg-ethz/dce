test_that("means are stable", {
  graph <- create_graph_from_dataframe(data.frame(
    from=c("A", "B"),
    to=c("B", "C")
  ), edge.weight=function() { 0 })

  X <- simulate_data(graph, n=1000, dist.mean=42, link.log.base=2)

  expect_equal(
    X %>% as.data.frame %>% summarise_all(list(mean)),
    data.frame(A=42, B=42, C=42),
    tolerance=0.01
  )
})

test_that("simple propagation works", {
  graph <- create_graph_from_dataframe(data.frame(
    from=c("A", "B"),
    to=c("B", "C")
  ), edge.weight=function() { 10 })

  X <- simulate_data(graph, n=1000, dist.mean=42, link.log.base=1.1)

  # TODO: this is wrong
  expect_equal(
    X %>% as.data.frame %>% summarise_all(list(mean)),
    data.frame(A=42, B=52, C=62),
    tolerance=0.01
  )
})
