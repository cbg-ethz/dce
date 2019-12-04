library(tidyverse)
library(furrr)
library(graph)

devtools::install("..", upgrade="never")
library(dce)


future::plan(multiprocess)
set.seed(42)


# parse commandline arguments
"
Benchmark DCE performance and runtime.

Usage:
  differential.R
  differential.R --variable NAME --values VALUES
" -> doc

arguments <- docopt::docopt(doc)


# parse parameters
node.num <- 100
wt.samples <- 200
mt.samples <- 200

varied.parameter <- "node.num"
parameter.list <- c(20, 50, 100)

if (!is.null(arguments$NAME)) {
  varied.parameter <- arguments$NAME
  parameter.list <- unlist(purrr::map(strsplit(arguments$VALUES, ",")[[1]], as.numeric))
}


# do benchmarking
replicate.count <- 100

df.bench <- furrr::future_pmap_dfr(
  list(parameter=rep(parameter.list, each=replicate.count)),
  function(parameter) {
    # handle parametrization
    switch(
      varied.parameter,
      node.num={ node.num <- parameter },
      wt.samples={ wt.samples <- parameter },
      mt.samples={ mt.samples <- parameter },
    )
    print(paste("Setting", varied.parameter, "to", parameter))


    # create graphs
    edge.prob <- runif(1, 0, 1)

    negweight.range <- c(-1, 0)
    posweight.range <- c(0, 1)

    wt.graph <- create_random_DAG(node.num, edge.prob, negweight.range, posweight.range)
    mt.graph <- resample_edge_weights(wt.graph, negweight.range, posweight.range)


    # compute graph features
    tmp <- as(wt.graph, "matrix")
    tmp[which(tmp != 0)] <- 1
    graph.density <- sum(tmp) / ((dim(tmp)[1] * (dim(tmp)[1] - 1)) / 2)


    # generate data
    wt.X <- simulate_data(wt.graph, n=wt.samples)
    mt.X <- simulate_data(mt.graph, n=mt.samples)


    # run models
    ground.truth <- list(dce=trueEffects(wt.graph) - trueEffects(mt.graph))

    time.tmp <- Sys.time()
    res.cor <- list(dce=cor(wt.X) - cor(mt.X))
    time.cor <- as.integer(difftime(Sys.time(), time.tmp, units="secs"))

    time.tmp <- Sys.time()
    res.basic <- compute_differential_causal_effects(
      wt.graph, wt.X,
      mt.graph, mt.X,
      method="basic"
    )
    time.basic <- as.integer(difftime(Sys.time(), time.tmp, units="secs"))

    time.tmp <- Sys.time()
    res.basic.bootstrap <- compute_differential_causal_effects(
      wt.graph, wt.X,
      mt.graph, mt.X,
      method="basic", bootstrap=TRUE
    )
    time.basic.bootstrap <- as.integer(difftime(Sys.time(), time.tmp, units="secs"))

    time.tmp <- Sys.time()
    res.full <- compute_differential_causal_effects(
      wt.graph, wt.X,
      mt.graph, mt.X,
      method="full"
    )
    time.full <- as.integer(difftime(Sys.time(), time.tmp, units="secs"))

    time.tmp <- Sys.time()
    res.full.bootstrap <- compute_differential_causal_effects(
      wt.graph, wt.X,
      mt.graph, mt.X,
      method="full", bootstrap=TRUE
    )
    time.full.bootstrap <- as.integer(difftime(Sys.time(), time.tmp, units="secs"))

    time.tmp <- Sys.time()
    tmp <- as.matrix(ground.truth$dce)
    tmp[which(as.matrix(ground.truth$dce) != 0)] = (
      runif(sum(as.matrix(ground.truth$dce) != 0), negweight.range[1], posweight.range[2]) -
      runif(sum(as.matrix(ground.truth$dce) != 0), negweight.range[1], posweight.range[2])
    )
    res.rand <- list(dce=tmp)
    time.rand <- as.integer(difftime(Sys.time(), time.tmp, units="secs"))

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

        bind_rows(list(
          cor=get_prediction_counts(df.res$truth, df.res$cor),
          basic=get_prediction_counts(df.res$truth, df.res$basic),
          basic.bootstrap=get_prediction_counts(df.res$truth, df.res$basic.bootstrap),
          full=get_prediction_counts(df.res$truth, df.res$full),
          full.bootstrap=get_prediction_counts(df.res$truth, df.res$full.bootstrap),
          rand=get_prediction_counts(df.res$truth, df.res$rand)
        ), .id="name") %>%
          column_to_rownames(var="name") %>%
          t %>%
          as.data.frame %>%
          rownames_to_column(var="type"),

        data.frame(
          cor=time.cor,
          basic=time.basic,
          basic.bootstrap=time.basic.bootstrap,
          full=time.full,
          full.bootstrap=time.full.bootstrap,
          rand=time.rand
        ) %>%
          mutate(type="runtime"),

        data.frame(
          cor=graph.density,
          basic=graph.density,
          basic.bootstrap=graph.density,
          full=graph.density,
          full.bootstrap=graph.density,
          rand=graph.density
        ) %>%
          mutate(type="graph.density"),
      ) %>%
      mutate(parameter=parameter)
  },
  .progress=TRUE
) %>%
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
  ggtitle(paste("Variable:", varied.parameter)) +
  ylab("performance") +
  theme_minimal() +
  theme(plot.title=element_text(hjust=0.5)) +
  ggsave("benchmark.pdf")

df.bench %>%
  dplyr::filter(type == "runtime") %>%
  gather("variable", "value", -parameter, -type) %>%
  mutate(value=lubridate::as.duration(value)) %>%
ggplot(aes(x=parameter, y=value, fill=variable)) +
  geom_boxplot() +
  scale_y_time() +
  ggtitle(paste("Variable:", varied.parameter)) +
  ylab("runtime") +
  theme_minimal() +
  theme(plot.title=element_text(hjust=0.5)) +
  ggsave("runtime.pdf")

df.bench %>%
  dplyr::filter(grepl("^graph.", type)) %>%
  select(full, type, parameter) %>%
ggplot(aes(x=parameter, y=full, fill=type)) +
  geom_boxplot() +
  ggtitle(paste("Variable:", varied.parameter)) +
  ylab("value") +
  theme_minimal() +
  theme(plot.title=element_text(hjust=0.5)) +
  ggsave("graph_features.pdf")
