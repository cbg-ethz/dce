library(tidyverse)
library(tidygraph)

library(pcalg)

# read data
df.expr <- read_csv("tcga_data/expression_matrix.csv")
df.classi <- read_csv("tcga_data/case_classifications.csv")
df.graph <- read_csv("kegg_data/hsa05219.edgelist.csv")

# remove suffix from Ensembl ids
df.expr.fltr <- df.expr %>%
  extract(gene, "gene.short") %>% 
  column_to_rownames("gene.short")


# compute causal effects
compute <- function (tissue.type) {
  print(tissue.type)
  
  # subset to relevant cases
  case.list <- df.classi %>%
    filter(tissue.definition == tissue.type) %>%
    pull(cases)
  node.list <- unique(c(df.graph$source, df.graph$sink))
  
  df.expr.sub <- t(df.expr.fltr[node.list, case.list])
  print(df.expr.sub %>% dim)
  
  # prepare data
  cov.mat <- cov(df.expr.sub)
  assertthat::are_equal(dim(cov.mat)[[1]], length(node.list))
  
  g <- tidygraph::as_tbl_graph(df.graph)
  ig <- igraph::as.igraph(g)
  gn <- igraph::igraph.to.graphNEL(ig)
  # pcalg::isValidGraph(gn)
  
  # estimate causal effects
  compute_causal_effect <- function (source, sink) {
    source.id <- match(source, node.list)
    sink.id <- match(sink, node.list)
    
    ida.res <- pcalg::ida(source.id, sink.id, cov.mat, gn, verbose=FALSE)
    return(ida.res)
  }
  
  df.graph %>% 
    rowwise() %>% 
    mutate(causal.effect=compute_causal_effect(source, sink))
}


compute("Solid Tissue Normal") %>%
  write_csv("results/graph.normal.edgelist.csv")

compute("Primary solid Tumor") %>%
  write_csv("results/graph.tumor.edgelist.csv")
