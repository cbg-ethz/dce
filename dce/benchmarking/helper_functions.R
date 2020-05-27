generate.random.graphs <- function(node.num, beta.magnitude, true.positives) {
  edge.prob <- runif(1, 0, 1)

  negweight.range <- c(-beta.magnitude, 0)
  posweight.range <- c(0, beta.magnitude)

  wt.graph <- create_random_DAG(node.num, edge.prob, negweight.range, posweight.range)
  while(length(wt.graph@edgeData@data) <= 1) {
    wt.graph <- create_random_DAG(node.num, edge.prob, negweight.range, posweight.range)
  }

  mt.graph <- resample_edge_weights(wt.graph, tp = true.positives,
                                    mineff = beta.magnitude,
                                    maxeff = 0.5,
                                    method = "exp")

  return(list(wt=wt.graph, mt=mt.graph))
}

sample.graph.from.kegg <- function(kegg.dag) {
  wt.graph <- kegg.dag[[sample(1:length(kegg.dag), 1)]]
  kegg.order <- order(apply(wt.graph, 1, sum) - apply(wt.graph, 2, sum), decreasing = TRUE)
  wt.graph <- wt.graph[kegg.order, kegg.order]
  wt.graph[lower.tri(wt.graph)] <- 0
  diag(wt.graph) <- 0
  wt.graph <- wt.graph[seq_len(min(node.num, nrow(wt.graph))), seq_len(min(node.num, nrow(wt.graph)))]
  wt.graph <- as(wt.graph, "graphNEL")
  wt.graph <- resample_edge_weights(wt.graph, negweight.range, posweight.range)

  mt.graph <- resample_edge_weights(wt.graph, negweight.range, posweight.range, tp = true.positives)

  return(list(wt=wt.graph, mt=mt.graph))
}

perturb.dag <- function(dag, perturb) {
  if (perturb > 0) {
    p.dag <- as(dag, "matrix")
    p.dag[which(p.dag != 0)] <- 1
    candidates <- intersect(which(p.dag == 0),
                            which(upper.tri(p.dag) == TRUE))
    turn <- sample(candidates, floor(perturb*length(candidates)))
    p.dag[turn] <- 1
    p.dag <- as(p.dag, "graphNEL")
  } else if (perturb < 0) {
    p.dag <- as(dag, "matrix")
    p.dag[which(p.dag != 0)] <- 1
    candidates <- which(p.dag != 0)
    sample.n <- floor(abs(perturb)*length(candidates))
    if (length(candidates) - sample.n <= 1) {
      sample.n <- length(candidates)
    }
    turn <- sample(candidates, sample.n)
    p.dag[turn] <- 0
    p.dag <- as(p.dag, "graphNEL")
  } else {
    p.dag <- dag
  }

  return(p.dag)
}

compute.prevalence <- function(wt.graph, mt.graph) {
  wt.a <- as(wt.graph, "matrix")
  mt.a <- as(mt.graph, "matrix")
  p <- sum(wt.a != mt.a)/sum(wt.a != 0)
  return(p)
}

compute.dce.stats <- function(wt.graph, mt.graph) {
  ground.truth <- abs(trueEffects(mt.graph) - trueEffects(wt.graph))
  ground.truth[which(ground.truth == 0)] <- NA
  res <- list(min = min(ground.truth, na.rm = TRUE),
              max = max(ground.truth, na.rm = TRUE),
              median = median(ground.truth, na.rm = TRUE),
              mean = mean(ground.truth, na.rm = TRUE))
  return(res)
}
