test_that("solver argument works", {
  set.seed(42)

  node_names <- c("A", "B")

  graph_wt <- matrix(c(0, 0, 1e-42, 0), 2, 2)
  rownames(graph_wt) <- colnames(graph_wt) <- node_names
  X_wt <- simulate_data(graph_wt)

  graph_mt <- matrix(c(0, 0, -1.5, 0), 2, 2)
  rownames(graph_mt) <- colnames(graph_mt) <- node_names
  X_mt <- simulate_data(graph_mt)

  custom_solver <- function(formula, data) {
    fit <- lm(formula = formula, data = data)
    summary(fit)$coefficients
  }

  res_1 <- dce::dce(graph_wt, X_wt, X_mt, solver = "lm")
  res_2 <- dce::dce(graph_wt, X_wt, X_mt, solver = lm)
  res_3 <- dce::dce(graph_wt, X_wt, X_mt, solver = custom_solver)
  res_4 <- dce::dce(graph_wt, X_wt, X_mt, solver = "glm2")
  res_5 <- dce::dce(graph_wt, X_wt, X_mt, solver = glm2::glm2)

  expect_equal(res_1, res_2)
  expect_equal(res_2, res_3)
  expect_equal(res_3, res_4, tolerance = 0.1)
  expect_equal(res_4, res_5)
})


test_that("positive beta can be recovered", {
  set.seed(42)

  node_names <- c("A", "B")

  graph_wt <- matrix(c(0, 0, 1e-42, 0), 2, 2)
  rownames(graph_wt) <- colnames(graph_wt) <- node_names
  X_wt <- simulate_data(graph_wt)

  graph_mt <- matrix(c(0, 0, 1.5, 0), 2, 2)
  rownames(graph_mt) <- colnames(graph_mt) <- node_names
  X_mt <- simulate_data(graph_mt)

  res <- dce::dce_nb(graph_wt, X_wt, X_mt)
  res

  expect_equal(as.vector(res$dce), c(NA, NA, 1.5, NA), tolerance = 0.1)

  expect_equal(rownames(res$dce), node_names)
  expect_equal(colnames(res$dce), node_names)
  expect_equal(rownames(res$dce.pvalue), node_names)
  expect_equal(colnames(res$dce.pvalue), node_names)
})


test_that("negative beta can be recovered", {
  set.seed(42)

  graph_wt <- matrix(c(0, 0, 1e-42, 0), 2, 2)
  X_wt <- simulate_data(graph_wt)

  graph_mt <- matrix(c(0, 0, -1.5, 0), 2, 2)
  X_mt <- simulate_data(graph_mt)

  res <- dce::dce_nb(graph_wt, X_wt, X_mt)
  res

  expect_equal(as.vector(res$dce), c(NA, NA, -1.5, NA), tolerance = 0.1)
})


test_that("igraph input works", {
  set.seed(42)

  graph <- igraph::make_tree(3, children = 3, mode = "out")
  graph_wt <- igraph::set.edge.attribute(graph, name = "weight", value = 1)
  graph_mt <- igraph::set.edge.attribute(graph, name = "weight", value = 2.4)

  X_wt <- simulate_data(graph_wt, n = 1000)
  X_mt <- simulate_data(graph_mt, n = 1000)

  res <- dce::dce_nb(graph, X_wt, X_mt)
  res

  expect_equal(
    as.vector(res$dce), c(NA, NA, NA, 1.4, NA, NA, 1.4, NA, NA),
    tolerance = 0.1
  )
})


test_that("graphNEL input works", {
  set.seed(42)

  graph_wt <- as(matrix(c(0, 0, 1e-42, 0), 2, 2), "graphNEL")
  X_wt <- simulate_data(graph_wt)

  graph_mt <- as(matrix(c(0, 0, 1.6, 0), 2, 2), "graphNEL")
  X_mt <- simulate_data(graph_mt)

  res <- dce::dce_nb(graph_wt, X_wt, X_mt)
  res

  expect_equal(as.vector(res$dce), c(NA, NA, 1.6, NA), tolerance = 0.1)
})


test_that("graph nodes and simulated data column mismatch throws error", {
  set.seed(42)

  graph <- dce::create_graph_from_dataframe(data.frame(
    from = c("A"),
    to = c("B")
  ))

  X_wt <- simulate_data(graph)
  X_mt <- simulate_data(graph)

  colnames(X_wt) <- c("Hinz", "Kunz") # introduce error

  expect_error(
    dce::dce_nb(graph, X_wt, X_mt),
    "Not all nodes have expression vector in WT data"
  )
})


test_that("adjustment sets work", {
  set.seed(42)

  graph <- dce::create_graph_from_dataframe(data.frame(
    from = c("A", "B"),
    to = c("B", "C")
  ))

  expect_equal(
    get_adjustment_set(
      as(graph, "matrix"),
      which(nodes(graph) == "B"), which(nodes(graph) == "C"),
      "parents"
    ),
    c("A")
  )
  expect_equal(
    get_adjustment_set(
      as(graph, "matrix"),
      which(nodes(graph) == "A"), which(nodes(graph) == "B"),
      "minimal"
    ),
    character(0)
  )
})


test_that("better solver can mitigate crash", {
  set.seed(42)

  # create problematic dataset
  beta <- 2
  theta <- 1
  mu <- 1000

  A <- rnbinom(1000, size = theta, mu = mu)

  beta.tmp <- beta * A
  B <- rnbinom(1000, size = theta, mu = beta.tmp - min(beta.tmp) + mu)

  # normal fit will fail
  expect_error(
    MASS::glm.nb(B ~ A, link = "identity", method = "glm.fit"),
    "no valid set of coefficients has been found: please supply starting values"
  )

  # custom fit will succeed
  fit <- MASS::glm.nb(B ~ A, link = "identity", method = "glm.dce.nb.fit")
  fit
})


test_that("plotting function does not crash", {
  set.seed(42)

  node_names <- c("A", "B", "C")

  graph_wt <- matrix(c(0, 0, 0, 1e-42, 0, 0, 1e-42, 0, 0), 3, 3)
  rownames(graph_wt) <- colnames(graph_wt) <- node_names
  X_wt <- simulate_data(graph_wt)

  graph_mt <- matrix(c(0, 0, 0, 1.5, 0, 0, -.7, 0, 0), 3, 3)
  rownames(graph_mt) <- colnames(graph_mt) <- node_names
  X_mt <- simulate_data(graph_mt)

  res <- dce::dce_nb(graph_wt, X_wt, X_mt)
  res

  plot(res)
})
