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
res <- dce::dce(
  igraph::induced_subgraph(graph, common.genes),
  X.wt, X.mt,
  solver = "lm",
  lib_size = FALSE, latent = 0
)

# save raw results
saveRDS(res, file = file.path(out.dir, glue::glue("dce_{appendix}.rds")))

# compute additional information
graph_sub <- igraph::induced_subgraph(graph, common.genes)

df_final <- res %>%
  as.data.frame %>%
  mutate(pathway_edge = melt(res$graph)$value) %>%
  filter(pathway_edge == 1) %>%
  select(-pathway_edge) %>%
  purrr::pmap_dfr(function(source, target, dce, dce_stderr, dce_pvalue) {
    dist <- igraph::distances(
      graph_sub,
      which(igraph::V(graph_sub)$name %in% strsplit(perturbed.gene, ",")[[1]]),
      which(igraph::V(graph_sub)$name %in% c(as.character(source), as.character(target))),
      mode = "all"
    ) %>%
      min

    data.frame(
      source = source, target = target,
      dce = dce, dce_stderr = dce_stderr, dce_pvalue = dce_pvalue,
      distance = dist
    )
  }) %>%
  arrange(dce_pvalue)

# save results as dataframe
df_final %>%
  write_csv(file.path(out.dir, glue::glue("dce_list_{appendix}.csv")))

# network plot
plot(res, labelsize = 1, highlighted_nodes = strsplit(perturbed.gene, ",")[[1]]) +
  ggtitle(glue::glue("{pathway}: {perturbed.gene}"))
ggsave(file.path(out.dir, glue::glue("network_{appendix}.pdf")), width = 20, height = 20)

# volcano plot
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

# p-value/network distance plot
df_final %>%
  filter(!is.infinite(distance)) %>%
ggplot(aes(x = -log10(dce_pvalue), y = distance, color = dce)) +
  geom_point() +
  xlab("-log10(DCE p-value)") +
  ylab("Graph distance perturbed-gene to DCE-edge") +
  ggtitle(glue::glue("{pathway}: {perturbed.gene}")) +
  theme_minimal()
ggsave(file.path(out.dir, glue::glue("network_distances_{appendix}.pdf")), width = 8, height = 6)
