library(tidyverse)

devtools::load_all("..")
source("../R/utils.R")


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


# create graph
node.num <- 30
edge.prob <- .2

negweight.range <- c(-1,-.1)
posweight.range <- c(.1, 1)

wt.graph <- randomDAG_2(node.num, edge.prob, negweight.range, posweight.range)
mt.graph <- newWeights(wt.graph, negweight.range, posweight.range, tp=1)

wt.adj <- t(as(wt.graph, "matrix"))
mt.adj <- t(as(mt.graph, "matrix"))


# do benchmarking
parameter.list <- c(10, 100, 1000)
repetition.num <- 20

tmp.list <- rep(parameter.list, each=repetition.num)
pb <- progress_estimated(length(tmp.list))
df.bench <- purrr::map_df(tmp.list, function(x) {
  pb$tick()$print()

  # generate data
  wt.X <- simulate(wt.graph, sample.num=x)
  mt.X <- simulate(mt.graph, sample.num=x)


  # run models
  ground.truth <- list(dce=trueEffects(wt.graph) - trueEffects(mt.graph))

  res.cor <- list(dce=cor(wt.X) - cor(mt.X))

  res.basic <- compute_differential_causal_effects(
    wt.graph, wt.X,
    mt.graph, mt.X,
    method="basic"
  )

  res.full <- compute_differential_causal_effects(
    wt.graph, wt.X,
    mt.graph, mt.X,
    method="full"
  )

  tmp <- as.matrix(ground.truth$dce)
  tmp[which(as.matrix(ground.truth$dce) != 0)] = (
    runif(sum(as.matrix(ground.truth$dce) != 0), negweight.range[1], posweight.range[2]) -
    runif(sum(as.matrix(ground.truth$dce) != 0), negweight.range[1], posweight.range[2])
  )
  res.rand <- list(dce=tmp)

  df.res <- data.frame(
    truth=as.vector(ground.truth$dce),
    cor=as.vector(res.cor$dce),
    basic=as.vector(res.basic$dce),
    full=as.vector(res.full$dce),
    rand=as.vector(res.rand$dce)
  )


  # compute performance
  df.perf <- as.data.frame(
    cor(df.res, method="spearman") #  %>% filter(truth != 0)
  )

  # return result
  df.perf$parameter <- as.factor(x)
  df.perf["truth", -1] # extract relations to "truth"
})

df.bench %>%
  head

df.bench %>%
  gather("variable", "value", -parameter) %>%
ggplot(aes(x=parameter, y=value, fill=variable)) +
  geom_boxplot() +
  ylim(-1, 1) +
  theme_minimal()

ggsave("benchmarking.pdf")
