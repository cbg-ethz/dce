library(tidyverse)

test_that("exponetial link function explodes mean", {
  set.seed(42)
  beta <- 0.1

  A <- rnbinom(100, size = 100, mu = 1000)
  B <- rnbinom(100, size = 100, mu = exp(A * beta))
  C <- rnbinom(100, size = 100, mu = exp(B * beta))

  data.frame(A=A, B=B, C=C) %>% summarize_all(list(mean))

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


test_that("CRISPR-like intervention leads to non-zero DCE", {
  set.seed(42)

  graph_wt <- dce::create_graph_from_dataframe(
    data.frame(
      from = c("A", "B"),
      to = c("B", "C")
    ),
    edge.weight = function() { 2 }
  )
  X_wt <- simulate_data(graph_wt, pop.size = 0)

  X_mt <- X_wt
  X_mt[, "B"] <- X_mt[, "B"] * 0.5

  res <- dce::dce_nb(graph_wt, X_wt, X_mt)
  res

  expect_lt(res$dce["A", "B"], -0.1)
})
