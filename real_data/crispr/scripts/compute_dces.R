devtools::load_all("../../../dce")

library(tidyverse)


# locate data
fname.graph <- snakemake@input$graph_file
fname.expr.wt <- snakemake@input$count_wt_file
fname.expr.mt <- snakemake@input$count_mt_file

out.dir <- snakemake@output$out_dir

perturbed.gene <- snakemake@wildcards$gene
appendix <- glue::glue("{snakemake@wildcards$pathway}_{perturbed.gene}")

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
#res <- dce::dce.nb(igraph::induced_subgraph(graph, common.genes), X.wt[, common.genes], X.mt[, common.genes])
res <- dce::dce.nb(igraph::induced_subgraph(graph, common.genes), X.wt, X.mt, lib.size = TRUE)

saveRDS(res, file = file.path(out.dir, glue::glue("dce_{appendix}.rds")))

# analyze results
plot(res, labelsize=1, highlighted.nodes=c(perturbed.gene))
ggsave(file.path(out.dir, glue::glue("network_{appendix}.pdf")), width=20, height=20)

res %>%
  as.data.frame %>%
  arrange(desc(abs(dce))) %>%
  drop_na %>%
  write_csv(file.path(out.dir, glue::glue("dce_list_{appendix}.csv")))
