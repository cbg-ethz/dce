library(tidyverse)
library(furrr)

devtools::install("..", upgrade="never")


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


# create graphs
graph.num <- 10

node.num <- 30
edge.prob <- .2

negweight.range <- c(-1,-.1)
posweight.range <- c(.1, 1)

graph.list <- purrr::imap(1:graph.num, function (x, i) {
  wt.graph <- randomDAG_2(node.num, edge.prob, negweight.range, posweight.range)
  mt.graph <- newWeights(wt.graph, negweight.range, posweight.range)

  list(graph.idx=i, wt.graph=wt.graph, mt.graph=mt.graph)
})

graph.list %>% head(1)


# compute graph features
graph.features <- purrr::map_df(graph.list, function (graph.pair) {
  tmp <- as(graph.pair$wt.graph, "matrix")
  tmp[which(tmp != 0)] <- 1

  data.frame(
    graph.idx=graph.pair$graph.idx,
    density=sum(tmp) / length(tmp)
  )
}) %>%
  write_csv("graph_features.csv")

graph.features %>% head(1)


# do benchmarking
parameter.list <- c(50, 100, 1000)
repetition.num <- 20

input.list <- purrr::cross_df(list(graph.pair=graph.list, parameter=parameter.list))

furrr::future_pmap_dfr(input.list, function(graph.pair, parameter) {
  purrr::map_df(seq_len(repetition.num), function(x) {
    # generate data
    wt.X <- simulate(graph.pair$wt.graph, sample.num=parameter)
    mt.X <- simulate(graph.pair$mt.graph, sample.num=parameter)


    # run models
    ground.truth <- list(dce=trueEffects(graph.pair$wt.graph) - trueEffects(graph.pair$mt.graph))

    time.tmp <- Sys.time()
    res.cor <- list(dce=cor(wt.X) - cor(mt.X))
    time.cor <- Sys.time() - time.tmp

    time.tmp <- Sys.time()
    res.basic <- compute_differential_causal_effects(
      graph.pair$wt.graph, wt.X,
      graph.pair$mt.graph, mt.X,
      method="basic"
    )
    time.basic <- Sys.time() - time.tmp

    time.tmp <- Sys.time()
    res.basic.bootstrap <- compute_differential_causal_effects(
      graph.pair$wt.graph, wt.X,
      graph.pair$mt.graph, mt.X,
      method="basic", bootstrap=TRUE
    )
    time.basic.bootstrap <- Sys.time() - time.tmp

    time.tmp <- Sys.time()
    res.full <- compute_differential_causal_effects(
      graph.pair$wt.graph, wt.X,
      graph.pair$mt.graph, mt.X,
      method="full"
    )
    time.full <- Sys.time() - time.tmp

    time.tmp <- Sys.time()
    res.full.bootstrap <- compute_differential_causal_effects(
      graph.pair$wt.graph, wt.X,
      graph.pair$mt.graph, mt.X,
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


    # compute performance
    df.perf <- as.data.frame(
      cor(df.res, method="spearman") #  %>% dplyr::filter(truth != 0)
    )

    # return result
    df.perf %>%
      rownames_to_column() %>%
      dplyr::filter(rowname == "truth") %>%
      select(-rowname, -truth) %>%
      mutate(type="performance") %>%
      bind_rows(
        data.frame(
          cor=time.cor,
          basic=time.basic,
          basic.bootstrap=time.basic.bootstrap,
          full=time.full,
          full.bootstrap=time.full.bootstrap,
          rand=time.rand
        ) %>%
          mutate(type="runtime")
      ) %>%
      mutate(parameter=parameter, graph.idx=graph.pair$graph.idx)
  })
}, .progress=TRUE) %>%
  write_csv("benchmark_results.csv")



df.bench$parameter <- fct_inseq(df.bench$parameter)
df.bench %>%
  head


# plotting
df.bench %>%
  dplyr::filter(type == "performance") %>%
  gather("variable", "value", -parameter, -type, -graph.idx) %>%
ggplot(aes(x=parameter, y=value, fill=variable)) +
  geom_boxplot() +
  ylim(-1, 1) +
  theme_minimal() +
  ggsave("benchmark.pdf")

df.bench %>%
  dplyr::filter(type == "runtime") %>%
  gather("variable", "value", -parameter, -type, -graph.idx) %>%
  mutate(value=lubridate::as.duration(value)) %>%
ggplot(aes(x=parameter, y=value, fill=variable)) +
  geom_boxplot() +
  scale_y_time() +
  theme_minimal() +
  ggsave("runtime.pdf")
