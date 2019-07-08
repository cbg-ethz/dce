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
df.between <- purrr::map_dfc(1:2, function (x) {
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
  res <- dce::compute_differential_causal_effects(graph, t(sub.normal), graph, t(sub.tumor))
  as.vector(res)
})

# within groups
get.within <- function (tissue.type) {
  purrr::map_dfc(1:2, function (x) {
    # sample subset
    cases.tumor1 <- df.classi %>%
      filter(tissue.definition == tissue.type) %>%
      sample_n(15) %>%
      pull(cases)
    cases.tumor2 <- df.classi %>%
      filter(tissue.definition == tissue.type) %>%
      sample_n(15) %>%
      pull(cases)

    sub.normal <- df.expr.fltr[node.list, cases.tumor1]
    sub.tumor <- df.expr.fltr[node.list, cases.tumor2]

    # compute
    res <- dce::compute_differential_causal_effects(graph, t(sub.normal), graph, t(sub.tumor))
    as.vector(res)
  })
}

df.within.normal <- get.within("Solid Tissue Normal")
df.within.tumor <- get.within("Primary solid Tumor")

# compare
df.comp <- data.frame(
  between=unlist(df.between),
  within.normal=unlist(df.within.normal),
  within.tumor=unlist(df.within.tumor)
) %>%
  gather("type", "dce") %>%
  write_csv(snakemake@output[["data_file"]])

ggplot(df.comp, aes(x=type, y=abs(dce))) +
  geom_violin() +
  scale_y_log10()
ggsave(snakemake@output[["plot_file"]])
