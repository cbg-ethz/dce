library(tidyverse)
library(magrittr)

devtools::load_all("../../../")


# read data
df.expr <- read_csv(snakemake@input$expression_fname) %>% column_to_rownames("gene")
df.classi <- read_csv(snakemake@input$classification_fname)
df.graph <- read_csv(snakemake@input$network_fname)

graph <- igraph::igraph.to.graphNEL(igraph::as.igraph(tidygraph::as_tbl_graph(df.graph)))
genes <- intersect(graph::nodes(graph), rownames(df.expr))


# compute stuff
tumor_stage_list <- c("stage i", "stage ii", "stage iii")

res <- purrr::map(tumor_stage_list, function(selected_tumor_stage) {
  # select control group
  barcodes.wt <- df.classi %>%
    dplyr::filter(definition == "Solid Tissue Normal") %>%
    pull(barcode)

  X.wt <- df.expr[, barcodes.wt] %>% t %>% as.data.frame


  # select experimental group
  barcodes.mt <- df.classi %>%
    dplyr::filter((definition == "Primary solid Tumor") & (tumor_stage == selected_tumor_stage)) %>%
    pull(barcode)

  X.mt <- df.expr[, barcodes.mt] %>% t %>% as.data.frame


  # annoying fix for nodes without data (TODO: improve this)
  undef.nodes <- setdiff(graph::nodes(graph), rownames(df.expr))
  X.wt[, undef.nodes] <- rnbinom(length(undef.nodes), size=1000, mu=100)
  X.mt[, undef.nodes] <- rnbinom(length(undef.nodes), size=1000, mu=100)


  # data statistics
  dim(X.wt)
  dim(X.mt)
  length(graph::nodes(graph))

  X.wt %>% summarize_all(list(mean))
  X.mt %>% summarize_all(list(mean))


  # compute DCEs
  dce::dce(graph, X.wt, X.mt, solver = "lm", test = "vcovHC", lib_size = TRUE)
}) %>%
  purrr::set_names(tumor_stage_list)


# store result
saveRDS(res, file = snakemake@output$dce_fname)

purrr::imap_dfr(res, function(x, name) {
  x %>%
    as.data.frame %>%
    drop_na %>%
    mutate(name = name)
}) %>%
  write_csv(snakemake@output$csv_fname)
