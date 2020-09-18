devtools::load_all("../../../dce")

library(tidyverse)


# locate data
fname.graph <- snakemake@input$graph_file
fname.expr.wt <- snakemake@input$count_wt_file
fname.expr.mt <- snakemake@input$count_mt_file

out.dir <- snakemake@output$out_dir
dir.create(out.dir, recursive = TRUE)

pathway <- snakemake@wildcards$pathway
perturbed.gene <- snakemake@wildcards$gene
appendix <- glue::glue("{pathway}_{perturbed.gene}")

# read data
X.wt <- read_csv(fname.expr.wt) %>%
  column_to_rownames("X1") %>%
  t
X.mt <- read_csv(fname.expr.mt) %>%
  column_to_rownames("X1") %>%
  t

graph <- igraph::graph.data.frame(read_csv(fname.graph))
common.genes <- intersect(igraph::vertex_attr(graph, "name"), colnames(X.wt))

# compute DCEs
#res <- dce::dce_nb(igraph::induced_subgraph(graph, common.genes), X.wt[, common.genes], X.mt[, common.genes])
res <- dce::dce_nb(igraph::induced_subgraph(graph, common.genes), X.wt, X.mt, lib_size = TRUE)

saveRDS(res, file = file.path(out.dir, glue::glue("dce_{appendix}.rds")))

# analyze results
res %>%
  as.data.frame %>%
  mutate(pathway_edge = melt(res$graph)$value) %>%
  filter(pathway_edge == 1) %>%
  select(-pathway_edge) %>%
  arrange(dce_pvalue) %>%
  write_csv(file.path(out.dir, glue::glue("dce_list_{appendix}.csv")))

plot(res, labelsize = 1, highlighted_nodes = strsplit(perturbed.gene, ",")[[1]])
ggsave(file.path(out.dir, glue::glue("network_{appendix}.pdf")), width = 20, height = 20)

df_volcano <- res %>%
  as.data.frame %>%
  drop_na %>%
  mutate(
    edge = paste0(source, "->", target),
    affected = (source %in% strsplit(perturbed.gene, ",")[[1]] | target %in% strsplit(perturbed.gene, ",")[[1]])
  )
if (dim(df_volcano)[[1]] > 0) {
  EnhancedVolcano::EnhancedVolcano(
    df_volcano,
    lab = df_volcano$edge, selectLab = filter(df_volcano, df_volcano$affected)$edge,
    x = "dce", y = "dce_pvalue",
    pCutoff = .05, FCcutoff = 1,
    drawConnectors = TRUE,
    title = glue::glue("{pathway}: {perturbed.gene}"), subtitle = NULL,
    xlab = bquote("DCE"), ylab = bquote(~-Log[10]~italic(pvalue)),
    legendLabels = c("NS", "DCE", "p-value", "p-value and DCE")
  )
  ggsave(file.path(out.dir, glue::glue("volcanoplot_{appendix}.pdf")), width = 10, height = 10)
}
