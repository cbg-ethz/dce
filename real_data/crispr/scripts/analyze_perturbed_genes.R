library(tidyverse)


# parameters
graph.files <- snakemake@params$all_pathways

perturbed.genes <- snakemake@params$perturbed_genes

out.dir <- snakemake@output$out_dir
dir.create(out.dir) # shouldn't snakemake do this automatically?

# handle multiple knockouts
perturbed.genes %<>%
  purrr::map(~ strsplit(., ",")[[1]]) %>%
  unlist %>%
  unique

# degree distribution
df.degree <- purrr::map_dfr(graph.files, function(fname) {
  graph <- igraph::graph.data.frame(read_csv(fname))
  pathway.name <- strsplit(basename(fname), "\\.")[[1]][[1]]

  genes <- intersect(igraph::vertex_attr(graph, "name"), perturbed.genes)
  if (length(genes) == 0) {
    return(data.frame())
  }

  degree.all <- igraph::degree(graph, mode = "all")[genes]
  degree.out <- igraph::degree(graph, mode = "out")[genes]
  degree.in <- igraph::degree(graph, mode = "in")[genes]

  data.frame(pathway = pathway.name, genes = genes, degree.all = degree.all, degree.out = degree.out, degree.in = degree.in)
}) %>%
  write_csv(file.path(out.dir, "perturbed_gene_statistics.csv"))

df.degree %>% arrange(desc(degree.in)) %>% head

ggplot(df.degree, aes(x=degree.in)) +
  geom_histogram() +
  theme_minimal()
ggsave(file.path(out.dir, "indegree_distribution.pdf"))
