library(tidyverse)
library(magrittr)


# read data
df.expr <- read_csv("../expr.csv.gz") %>% column_to_rownames("gene")
df.classi <- read_csv("../info.csv.gz")
df.graph <- read_csv("../netw.edgelist.csv.gz")

df.expr <- read_csv(snakemake@input$expression_fname) %>% column_to_rownames("gene")
df.classi <- read_csv(snakemake@input$classification_fname)
df.graph <- read_csv(snakemake@input$network_fname)

graph <- igraph::igraph.to.graphNEL(igraph::as.igraph(tidygraph::as_tbl_graph(df.graph)))


# group data
barcodes.wt <- df.classi %>%
  dplyr::filter(definition == "Solid Tissue Normal") %>%
  pull(barcode)

barcodes.mt <- df.classi %>%
  dplyr::filter(definition == "Primary solid Tumor") %>%
  pull(barcode)

genes <- intersect(nodes(graph), rownames(df.expr))
X.wt <- df.expr[genes, barcodes.wt] %>% t %>% as.data.frame
X.mt <- df.expr[genes, barcodes.mt] %>% t %>% as.data.frame


# annoying fix for nodes without data (TODO: improve this)
undef.nodes <- setdiff(nodes(graph), rownames(df.expr))
X.wt[, undef.nodes] <- rnbinom(length(undef.nodes), size=1000, mu=100)
X.mt[, undef.nodes] <- rnbinom(length(undef.nodes), size=1000, mu=100)


# uh oh (TODO: remove this)
m <- as(graph, "matrix")
ord <- gRbase::topo_sort(m, index=TRUE)
graph <- as(m[ord, ord], "graphNEL")

X.wt <- dce::simulate_data(graph) %>% as.data.frame
X.mt <- dce::simulate_data(graph) %>% as.data.frame
# uh oh end


# data statistics
dim(X.wt)
dim(X.mt)
length(nodes(graph))

X.wt %>% summarize_all(list(mean))
X.mt %>% summarize_all(list(mean))


# compute DCEs
res <- dce::compute_differential_causal_effects(
  graph, X.wt, 
  graph, X.mt
)


# store result
save(res, file=snakemake@output$dce_fname)