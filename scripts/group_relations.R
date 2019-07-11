library(tidyverse)
library(ggplot2)

library(dce)


# read data
df.expr <- read_csv(snakemake@input[["expression_file"]])
df.classi <- read_csv(snakemake@input[["classification_file"]])
df.graph <- read_csv(snakemake@input[["network_file"]])

g <- tidygraph::as_tbl_graph(df.graph)
ig <- igraph::as.igraph(g)
graph <- igraph::igraph.to.graphNEL(ig)

node.list <- unique(c(df.graph$source, df.graph$sink))

# remove suffix from Ensembl ids
df.expr.fltr <- df.expr %>%
  extract(gene, "gene.short") %>%
  column_to_rownames("gene.short")

# basic test
sub1 <- df.expr.fltr[node.list, sample(colnames(df.expr.fltr), 30)]
sub2 <- df.expr.fltr[node.list, sample(colnames(df.expr.fltr), 30)]

# debugonce(compute_causal_effects)
dce::compute_differential_causal_effects(graph, t(sub1), graph, t(sub2))

# between groups
df.between <- purrr::map_df(1:2, function (x) {
  # sample subset
  cases.normal <- df.classi %>%
    filter(tissue.definition == "Solid Tissue Normal") %>%
    sample_n(15) %>%
    pull(cases)
  cases.tumor <- df.classi %>%
    filter(tissue.definition == "Primary solid Tumor") %>%
    sample_n(15) %>%
    pull(cases)

  sub.normal <- df.expr.fltr[node.list, cases.normal]
  sub.tumor <- df.expr.fltr[node.list, cases.tumor]

  # compute
  ## DCE
  res <- dce::compute_differential_causal_effects(graph, t(sub.normal), graph, t(sub.tumor))

  ## DGE (TODO: use all genes for analysis?)
  sub.agg <- cbind(sub.normal, sub.tumor)
  group <- c(rep("N", length(cases.normal)), rep("T", length(cases.tumor)))
  y <- edgeR::DGEList(counts=sub.agg, group=group)
  y <- edgeR::calcNormFactors(y)
  design <- model.matrix(~group)
  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmQLFit(y, design)
  qlf <- edgeR::glmQLFTest(fit)#, coef=2)
  tt <- edgeR::topTags(qlf, n=NULL)$table

  # return data
  res.vec <- as.vector(res)
  data.frame(
    iteration=x,
    comparison.type="between",
    value.type=c(rep("dce", length(res.vec)), rep("dge", length(tt$logFC))),
    value=c(res.vec, tt$logFC)
  )
})

# within groups
get.within <- function (tissue.type) {
  purrr::map_df(1:2, function (x) {
    # sample subset
    cases.class1 <- df.classi %>%
      filter(tissue.definition == tissue.type) %>%
      sample_n(15) %>%
      pull(cases)
    cases.class2 <- df.classi %>%
      filter(tissue.definition == tissue.type) %>%
      sample_n(15) %>%
      pull(cases)

    sub.class1 <- df.expr.fltr[node.list, cases.class1]
    sub.class2 <- df.expr.fltr[node.list, cases.class2]

    # compute
    ## DCE
    res <- dce::compute_differential_causal_effects(graph, t(sub.class1), graph, t(sub.class2))

    ## DGE (TODO: use all genes for analysis?)
    sub.agg <- cbind(sub.class1, sub.class2)
    group <- c(rep("C1", length(sub.class1)), rep("C2", length(sub.class2)))
    y <- edgeR::DGEList(counts=sub.agg, group=group)
    y <- edgeR::calcNormFactors(y)
    design <- model.matrix(~group)
    y <- edgeR::estimateDisp(y, design)
    fit <- edgeR::glmQLFit(y, design)
    qlf <- edgeR::glmQLFTest(fit)#, coef=2)
    tt <- edgeR::topTags(qlf, n=NULL)$table

    # return data
    res.vec <- as.vector(res)
    data.frame(
      iteration=x,
      comparison.type=paste("within", tissue.type),
      value.type=c(rep("dce", length(res.vec)), rep("dge", length(tt$logFC))),
      value=c(res.vec, tt$logFC)
    )
  })
}

df.within.normal <- get.within("Solid Tissue Normal")
df.within.tumor <- get.within("Primary solid Tumor")

# compare
df.comp <- rbind(df.between, df.within.normal, df.within.tumor) %>%
  write_csv(snakemake@output[["data_file"]])

ggplot(df.comp, aes(x=comparison.type, y=abs(value))) +
  geom_violin() +
  scale_y_log10() +
  facet_grid(value.type ~ .)
ggsave(snakemake@output[["plot_file"]])
