library(tidyverse)
library(furrr)
library(graph)

devtools::install("..", upgrade="never")
library(dce)


future::plan(multiprocess)
set.seed(42)


# helper functions
simulate <- function(graph, noise.sd=1, sample.num=100) {
  p <- length(nodes(graph))
  adj <- t(as(graph, "matrix"))

  X <- matrix(0, nrow=sample.num, ncol=p)
  colnames(X) <- nodes(graph)
  X[, 1] <- 1 + rnorm(sample.num, mean=0, sd=noise.sd) # activate root node

  for (j in 2:p) {
    ij <- 1:(j - 1)
    X[, j] <- X[, j] + X[, ij, drop=FALSE] %*% adj[j, ij] + rnorm(sample.num, mean=0, sd=noise.sd)
  }

  X
}


# do benchmarking
node.num <- 30

parameter.list <- c(50, 100, 1000)

df.bench <- furrr::future_pmap_dfr(list(parameter=parameter.list), function(parameter) {
  # create graphs
  edge.prob <- runif(1, 0, 1)

  negweight.range <- c(-1, 0)
  posweight.range <- c(0, 1)

  wt.graph <- create_random_DAG(node.num, edge.prob, negweight.range, posweight.range)
  mt.graph <- resample_edge_weights(wt.graph, negweight.range, posweight.range)


  # generate data
  wt.X <- simulate(wt.graph, sample.num=parameter)
  mt.X <- simulate(mt.graph, sample.num=parameter)


  # run models
  ground.truth <- list(dce=trueEffects(wt.graph) - trueEffects(mt.graph))

  time.tmp <- Sys.time()
  res.cor <- list(dce=cor(wt.X) - cor(mt.X))
  time.cor <- Sys.time() - time.tmp

  time.tmp <- Sys.time()
  res.basic <- compute_differential_causal_effects(
    wt.graph, wt.X,
    mt.graph, mt.X,
    method="basic"
  )
  time.basic <- Sys.time() - time.tmp

  time.tmp <- Sys.time()
  res.basic.bootstrap <- compute_differential_causal_effects(
    wt.graph, wt.X,
    mt.graph, mt.X,
    method="basic", bootstrap=TRUE
  )
  time.basic.bootstrap <- Sys.time() - time.tmp

  time.tmp <- Sys.time()
  res.full <- compute_differential_causal_effects(
    wt.graph, wt.X,
    mt.graph, mt.X,
    method="full"
  )
  time.full <- Sys.time() - time.tmp

  time.tmp <- Sys.time()
  res.full.bootstrap <- compute_differential_causal_effects(
    wt.graph, wt.X,
    mt.graph, mt.X,
    method="full", bootstrap=TRUE
  )
  time.full.bootstrap <- Sys.time() - time.tmp

  time.tmp <- Sys.time()
  tmp <- as.matrix(ground.truth$dce)
  tmp[which(as.matrix(ground.truth$dce) != 0)] = (
    runif(sum(as.matrix(ground.truth$dce) != 0), negweight.range[1], posweight.range[2]) -
    runif(sum(as.matrix(ground.truth$dce) != 0), negweight.range[1], posweight.range[2])
  )
  res.rand <- list(dce=tmp)
  time.rand <- Sys.time() - time.tmp

  df.res <- data.frame(
    truth=as.vector(ground.truth$dce),
    cor=as.vector(res.cor$dce),
    basic=as.vector(res.basic$dce),
    basic.bootstrap=as.vector(res.basic.bootstrap$dce),
    full=as.vector(res.full$dce),
    full.bootstrap=as.vector(res.full.bootstrap$dce),
    rand=as.vector(res.rand$dce)
  )


  # return performance computation
  data.frame() %>%
    bind_rows(
      as.data.frame(
        cor(df.res, method="spearman")
      ) %>%
        rownames_to_column() %>%
        dplyr::filter(rowname == "truth") %>%
        select(-rowname, -truth) %>%
        mutate(type="correlation"),

      data.frame(
        cor=time.cor,
        basic=time.basic,
        basic.bootstrap=time.basic.bootstrap,
        full=time.full,
        full.bootstrap=time.full.bootstrap,
        rand=time.rand
      ) %>%
        mutate(type="runtime"),

    ) %>%
    mutate(parameter=parameter)
}, .progress=TRUE) %>%
  write_csv("benchmark_results.csv")


df.bench$parameter %<>% as.factor %>% fct_inseq
df.bench %>%
  head


# plotting
df.bench %>%
  dplyr::filter(type == "correlation") %>%
  gather("variable", "value", -parameter, -type) %>%
ggplot(aes(x=parameter, y=value, fill=variable)) +
  geom_boxplot() +
  ylim(-1, 1) +
  theme_minimal() +
  ggsave("benchmark.pdf")

df.bench %>%
  dplyr::filter(type == "runtime") %>%
  gather("variable", "value", -parameter, -type) %>%
  mutate(value=lubridate::as.duration(value)) %>%
ggplot(aes(x=parameter, y=value, fill=variable)) +
  geom_boxplot() +
  scale_y_time() +
  theme_minimal() +
  ggsave("runtime.pdf")
