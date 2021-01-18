test_that("exponential link function explodes mean", {
  set.seed(42)
  beta <- 0.2

  A <- rnbinom(100, size = 100, mu = 1000)
  B <- rnbinom(100, size = 100, mu = exp(A * beta))
  C <- rnbinom(100, size = 100, mu = exp(B * beta))

  data.frame(A = A, B = B, C = C) %>%
    summarize_all(list(mean))

  expect_true(all(is.na(C)))
})


test_that("library size correction is useful", {
  set.seed(42)
  N <- 100
  beta <- 1.6

  # generate data
  linkfun <- dce::negative.binomial.special()$linkfun
  A <- rnbinom(N, size = 100, mu = 1000)
  B <- rnbinom(N, size = 100, mu = linkfun(beta * A))

  # make library size correction necessary
  lib.factor <- sample(c(1, 2, 4, 8), 100, replace = TRUE)

  A.s <- A * lib.factor
  B.s <- B * lib.factor

  bg <- matrix(rnbinom(1000 * N, size = 100, mu = 1000), 1000, N)
  bg.s <- t(t(bg) * lib.factor)
  lib.size <- apply(bg.s, 2, sum)

  lib.size <- round(lib.size/(10^min(round(log10(lib.size)))))

  # data overview
  plot(data.frame(A=A, B=B, A.s=A.s, B.s=B.s))

  # helper function
  glm.fun <- function(...) {
    dce::glm.nb.rob(..., method = "glm.dce.nb.fit", link = "identity")
  }

  # fit models
  glm.fun(B ~ A) # works fine (no library size effect)
  glm.fun(B.s ~ A.s) # yields wrong estimate
  glm.fun(B.s ~ A.s + factor(lib.factor)) # works fine (but `lib.factor` is ground truth)
  glm.fun(B.s ~ A.s + factor(lib.size)) # uses realistic library size correction and yields reasonable estimate
})


test_that("transitive edges can affect total pathway enrichment", {
  set.seed(42)

  node_names <- c("A", "B")

  graph_wt <- matrix(c(0, 0, 1e-42, 0), 2, 2)
  rownames(graph_wt) <- colnames(graph_wt) <- node_names
  X_wt <- simulate_data(graph_wt)

  graph_mt <- matrix(c(0, 0, .2, 0), 2, 2)
  rownames(graph_mt) <- colnames(graph_mt) <- node_names
  X_mt <- simulate_data(graph_mt)

  res <- dce::dce_nb(graph_wt, X_wt, X_mt)
  res

  # TODO: finish this
})


test_that("CRISPR-like intervention leads to non-zero DCE", {
  set.seed(42)

  # positive beta leads to negative DCE
  beta <- 2
  graph_wt <- dce::create_graph_from_dataframe(
    data.frame(
      from = c("A", "B"),
      to = c("B", "C")
    ),
    edge_weight = function() { beta }
  )
  X_wt <- simulate_data(graph_wt, n = 100)

  X_mt <- X_wt
  X_mt[, "B"] <- X_mt[, "B"] * 0.5
  X_mt[, "C"] <- rnbinom(100, size = Inf, mu = dce::negative.binomial.special()$linkfun(X_mt[, "B"] * beta))

  res <- dce::dce_nb(graph_wt, X_wt, X_mt)
  res %>%
    as.data.frame

  plot(res)

  expect_lt(res$dce["A", "B"], 0)


  # negative beta leads to positive DCE
  beta <- -2
  graph_wt <- dce::create_graph_from_dataframe(
    data.frame(
      from = c("A", "B"),
      to = c("B", "C")
    ),
    edge_weight = function() { beta }
  )
  X_wt <- simulate_data(graph_wt, n = 100)

  X_mt <- X_wt
  X_mt[, "B"] <- X_mt[, "B"] * 0.5
  X_mt[, "C"] <- rnbinom(100, size = Inf, mu = dce::negative.binomial.special()$linkfun(X_mt[, "B"] * beta))

  res <- dce::dce_nb(graph_wt, X_wt, X_mt)
  res %>%
    as.data.frame

  plot(res)

  expect_gt(res$dce["A", "B"], 0)
})
