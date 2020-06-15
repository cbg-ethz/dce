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
