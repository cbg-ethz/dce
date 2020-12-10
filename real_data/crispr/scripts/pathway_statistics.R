library(tidyverse)


# parameters
fname.expr.wt <- snakemake@input$count_wt_file
graph.files <- snakemake@input$graph_files

perturbed.genes <- snakemake@params$perturbed_genes

out.dir <- snakemake@output$out_dir
dir.create(out.dir, recursive = TRUE)

# handle multiple knockouts
perturbed.genes %<>%
  purrr::map(~ strsplit(., ",")[[1]]) %>%
  unlist %>%
  unique

# degree distribution
df.degree <- purrr::map_dfr(graph.files, function(fname) {
  graph <- igraph::graph.data.frame(read_csv(
    fname, col_types = cols(
      source = col_character(),
      sink = col_character()
    )
  ))
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
  arrange(desc(degree.all)) %>%
  write_csv(file.path(out.dir, "perturbed_gene_pathway_degrees.csv"))

df.degree %>% arrange(desc(degree.in)) %>% head

ggplot(df.degree, aes(x = degree.in)) +
  geom_histogram() +
  theme_minimal()
ggsave(file.path(out.dir, "indegree_distribution.pdf"))

# count distribution
df.expr.wt <- read_csv(fname.expr.wt) %>% column_to_rownames("X1") %>% as.matrix

df.counts <- purrr::map_dfr(graph.files, function(fname) {
  graph <- igraph::graph.data.frame(read_csv(
    fname, col_types = cols(
      source = col_character(),
      sink = col_character()
    )
  ))
  pathway.name <- strsplit(basename(fname), "\\.")[[1]][[1]]

  pathway_nodes <- igraph::vertex_attr(graph, "name")
  common_genes <- intersect(pathway_nodes, rownames(df.expr.wt))

  df.sub <- df.expr.wt[common_genes, ]
  if (is.null(dim(df.sub))) {
    # only a single gene was selected
    stopifnot(length(common_genes) == 1)

    gene_medians <- median(df.sub)
  } else {
    gene_medians <- matrixStats::rowMedians(df.sub)
  }

  good_counts <- sum(gene_medians > 1)
  good_count_fraction <- good_counts / length(gene_medians)

  data.frame(
    pathway = pathway.name,
    pathway_size = length(pathway_nodes),
    pathway_countmatrix_overlap = length(common_genes),
    good_counts = good_counts,
    good_count_fraction = good_count_fraction
  )
}) %>%
  arrange(desc(good_count_fraction)) %>%
  write_csv(file.path(out.dir, "pathway_count_statistics.csv"))

ggplot(df.counts, aes(x = good_count_fraction)) +
  geom_histogram() +
  xlab("Fraction of pathway nodes with median counts > 1") +
  theme_minimal()
ggsave(file.path(out.dir, "good_count_distribution.pdf"))
