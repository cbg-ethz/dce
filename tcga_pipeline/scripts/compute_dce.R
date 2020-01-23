library(tidyverse)
library(magrittr)


# read data
df.expr <- read_csv(snakemake@input$expression_fname) %>% column_to_rownames("gene")
df.classi <- read_csv(snakemake@input$classification_fname)
df.graph <- read_csv(snakemake@input$network_fname)

graph <- igraph::igraph.to.graphNEL(igraph::as.igraph(tidygraph::as_tbl_graph(df.graph)))
genes <- intersect(graph::nodes(graph), rownames(df.expr))


# compute stuff
tumor_stage_list <- c("stage i", "stage ii", "stage iii")

res <- purrr::map(tumor_stage_list, function (selected_tumor_stage) {
  # select control group
  barcodes.wt <- df.classi %>%
    dplyr::filter(definition == "Solid Tissue Normal") %>%
    pull(barcode)

  X.wt <- df.expr[genes, barcodes.wt] %>% t %>% as.data.frame


  # select experimental group
  barcodes.mt <- df.classi %>%
    dplyr::filter((definition == "Primary solid Tumor") & (tumor_stage == selected_tumor_stage)) %>%
    pull(barcode)

  X.mt <- df.expr[genes, barcodes.mt] %>% t %>% as.data.frame


  # annoying fix for nodes without data (TODO: improve this)
  undef.nodes <- setdiff(graph::nodes(graph), rownames(df.expr))
  X.wt[, undef.nodes] <- rnbinom(length(undef.nodes), size=1000, mu=100)
  X.mt[, undef.nodes] <- rnbinom(length(undef.nodes), size=1000, mu=100)


  # uh oh (TODO: remove this)
  m <- as(graph, "matrix")
  ord <- gRbase::topo_sort(m, index=TRUE)
  graph <- as(m[ord, ord], "graphNEL")

  X.wt <- dce::simulate_data(graph) %>% as.data.frame

  resample.range <- list(`stage i`=c(0.5, 1.5), `stage ii`=c(1.5, 2), `stage iii`=c(2, 3))[[selected_tumor_stage]]
  X.mt <- dce::simulate_data(dce::resample_edge_weights(graph, lB=resample.range, uB=resample.range)) %>% as.data.frame
  # uh oh end


  # data statistics
  dim(X.wt)
  dim(X.mt)
  length(graph::nodes(graph))

  X.wt %>% summarize_all(list(mean))
  X.mt %>% summarize_all(list(mean))


  # compute DCEs
  dce::compute_differential_causal_effects(
    graph, X.wt,
    graph, X.mt
  )
}) %>%
  purrr::set_names(tumor_stage_list)


# store result
save(res, file=snakemake@output$dce_fname)
