library(tidyverse)


# parameters
fname_expr_wt <- snakemake@input$count_wt_file
graph_files <- snakemake@input$graph_files

perturbed_genes <- snakemake@params$perturbed_genes

out_dir <- snakemake@output$out_dir
dir.create(out_dir, recursive = TRUE)


# handle multiple knockouts
perturbed_genes %<>%
  purrr::map(~ strsplit(., ",")[[1]]) %>%
  unlist %>%
  unique


# read auxiliary data
df_expr_wt <- read_csv(fname_expr_wt) %>% column_to_rownames("X1") %>% as.matrix


# compute per-pathway statistics
df_stats <- purrr::map_dfr(graph_files, function(fname) {
  # read pathway
  graph <- igraph::graph.data.frame(read_csv(
    fname, col_types = cols(
      source = col_character(),
      sink = col_character()
    )
  ))
  pathway_name <- strsplit(basename(fname), "\\.")[[1]][[1]]

  # pathway graph properties
  mean_degree <- mean(igraph::degree(graph, mode = "all"))
  mean_degree_out <- mean(igraph::degree(graph, mode = "out"))
  mean_degree_in <- mean(igraph::degree(graph, mode = "in"))

  # determine expression coverage along pathway
  pathway_nodes <- igraph::vertex_attr(graph, "name")
  common_genes <- intersect(pathway_nodes, rownames(df_expr_wt))

  df_sub <- df_expr_wt[common_genes, ]
  if (is.null(dim(df_sub))) {
    # only a single gene was selected
    stopifnot(length(common_genes) == 1)

    gene_medians <- median(df_sub)
  } else {
    gene_medians <- matrixStats::rowMedians(df_sub)
  }

  meanofmedian_gene_expression <- mean(gene_medians)

  data.frame(
    pathway = pathway_name,

    pathway_node_size = igraph::vcount(graph),
    pathway_edge_size = igraph::ecount(graph),
    pathway_density = igraph::ecount(graph) / (igraph::vcount(graph) * (igraph::vcount(graph) - 1)),
    mean_degree = mean_degree,
    mean_degree_out = mean_degree_out,
    mean_degree_in = mean_degree_in,

    pathway_countmatrix_overlap = length(common_genes),

    meanofmedian_gene_expression = meanofmedian_gene_expression
  )
}) %>%
  write_csv(file.path(out_dir, "pathway_statistics.csv"))


# degree distribution of perturbed genes
df_degree <- purrr::map_dfr(graph_files, function(fname) {
  graph <- igraph::graph.data.frame(read_csv(
    fname, col_types = cols(
      source = col_character(),
      sink = col_character()
    )
  ))
  pathway_name <- strsplit(basename(fname), "\\.")[[1]][[1]]

  genes <- intersect(igraph::vertex_attr(graph, "name"), perturbed_genes)
  if (length(genes) == 0) {
    return(data.frame())
  }

  degree_all <- igraph::degree(graph, mode = "all")[genes]
  degree_out <- igraph::degree(graph, mode = "out")[genes]
  degree_in <- igraph::degree(graph, mode = "in")[genes]

  data.frame(
    pathway = pathway_name,
    genes = genes,

    degree_all = degree_all,
    degree_out = degree_out,
    degree_in = degree_in
  )
}) %>%
  arrange(desc(degree_all)) %>%
  write_csv(file.path(out_dir, "perturbed_gene_pathway_degrees.csv"))

df_degree %>% arrange(desc(degree_in)) %>% head

ggplot(df_degree, aes(x = degree_in)) +
  geom_histogram() +
  theme_minimal()
ggsave(file.path(out_dir, "indegree_distribution.pdf"))


# pathway statistic plots
df_stats %>%
  select_if(is.numeric) %>%
  gather %>%
ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~ key, scales = "free") +
  theme_minimal()
ggsave(file.path(out_dir, "pathway_statistics.pdf"))
