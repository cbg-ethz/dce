library(tidyverse)

devtools::load_all("../../../")

# locate data
fname.graph <- snakemake@input$graph_file
fname.expr.wt <- snakemake@input$count_wt_file
fname.expr.mt <- snakemake@input$count_mt_file

out.dir <- snakemake@output$out_dir
dir.create(out.dir, recursive = TRUE)

pathway <- snakemake@wildcards$pathway
perturbed.gene <- snakemake@wildcards$gene
appendix <- glue::glue("{pathway}_{perturbed.gene}")

params <- snakemake@params
if (params$computation$deconfounding == "is_true") {
  params$computation$deconfounding <- TRUE
} else if (params$computation$deconfounding == "is_false") {
  params$computation$deconfounding <- FALSE
} else {
  conv <- as.numeric(params$computation$deconfounding)
  if (!is.na(conv)) {
    # parameter is a number (otherwise string of method name)
    params$computation$deconfounding <- conv
  }
}
print(params)

# read data
X.wt <- read_csv(fname.expr.wt) %>%
  column_to_rownames_wrap("...1") %>%
  t
X.mt <- read_csv(fname.expr.mt) %>%
  column_to_rownames_wrap("...1") %>%
  t

graph <- igraph::graph.data.frame(read_csv(fname.graph))
common.genes <- intersect(intersect(igraph::vertex_attr(graph, "name"), colnames(X.wt)), colnames(X.mt))

if (length(common.genes) > 0) {
  
  # compute DCEs
  res <- dce::dce(
    igraph::induced_subgraph(graph, common.genes),
    X.wt[,common.genes], X.mt[,common.genes],
    solver = "lm",
    test = "vcovHC",
    deconfounding = params$computation$deconfounding
  )
  
  # run competing models
  res_cor <- list(
    dce = cor(X.mt[, common.genes]) - cor(X.wt[, common.genes])
  )
  res_cor$dce_pvalue <- dce::permutation_test(
    X.wt[, common.genes], X.mt[, common.genes],
    fun = cor
  )
  
  res_pcor <- list(
    dce = dce::pcor(X.mt[, common.genes]) - dce::pcor(X.wt[, common.genes])
  )
  res_pcor$dce_pvalue <- dce::permutation_test(
    X.wt[, common.genes], X.mt[, common.genes],
    fun = dce::pcor, iter = 10
  )
  
  # clean p-values
  min_val <- min(res_cor$dce_pvalue[res_cor$dce_pvalue != 0], na.rm = TRUE)
  res_cor$dce_pvalue[res_cor$dce_pvalue == 0] <- min_val
  
  min_val <- min(res_pcor$dce_pvalue[res_pcor$dce_pvalue != 0], na.rm = TRUE)
  res_pcor$dce_pvalue[res_pcor$dce_pvalue == 0] <- min_val
  
  
  # plot method comparison
  p <- cowplot::plot_grid(
    plot_network(
      res$graph, value_matrix = res$dce,
      legend_title = "DCE",
      labelsize = 1, highlighted_nodes = strsplit(perturbed.gene, ",")[[1]]
    ),
    plot_network(
      res$graph, value_matrix = -log10(res$dce_pvalue),
      legend_title = "DCE p-value (-log)",
      edgescale_limits = c(0, max(-log10(res$dce_pvalue), na.rm = TRUE)),
      labelsize = 1, highlighted_nodes = strsplit(perturbed.gene, ",")[[1]]
    ),
    plot_network(
      res$graph, value_matrix = res_cor$dce,
      legend_title = "Corr",
      labelsize = 1, highlighted_nodes = strsplit(perturbed.gene, ",")[[1]]
    ),
    plot_network(
      res$graph, value_matrix = -log10(res_cor$dce_pvalue),
      legend_title = "Corr p-value (-log)",
      edgescale_limits = c(0, max(-log10(res_cor$dce_pvalue), na.rm = TRUE)),
      labelsize = 1, highlighted_nodes = strsplit(perturbed.gene, ",")[[1]]
    ),
    plot_network(
      res$graph, value_matrix = res_pcor$dce,
      legend_title = "P-Corr",
      labelsize = 1, highlighted_nodes = strsplit(perturbed.gene, ",")[[1]]
    ),
    plot_network(
      res$graph, value_matrix = -log10(res_pcor$dce_pvalue),
      legend_title = "P-Corr p-value (-log)",
      edgescale_limits = c(0, max(-log10(res_pcor$dce_pvalue), na.rm = TRUE)),
      labelsize = 1, highlighted_nodes = strsplit(perturbed.gene, ",")[[1]]
    ),
    nrow = 3, ncol = 2,
    labels = c(
      "DCE", "DCE p-value",
      "Corr", "Corr p-value",
      "P-Corr", "P-Corr p-value"
    )
  )
  p
  cowplot::save_plot(
    filename = file.path(out.dir, glue::glue("method_comparison_{appendix}.pdf")),
    plot = p,
    nrow = 3, ncol = 2,
    base_height = 8, base_asp = 1, limitsize = FALSE
  )
  
  # save raw results
  saveRDS(res, file = file.path(out.dir, glue::glue("dce_{appendix}.rds")))
  
  # compute additional information
  graph_sub <- igraph::induced_subgraph(graph, common.genes)
  
  df_final <- res %>%
    as.data.frame %>%
    
    # add competing methods
    mutate(cor = melt(res_cor$dce)$value) %>%
    mutate(cor_pvalue = melt(res_cor$dce_pvalue)$value) %>%
    mutate(pcor = melt(res_pcor$dce)$value) %>%
    mutate(pcor_pvalue = melt(res_pcor$dce_pvalue)$value) %>%
    
    # retain only values on pathway edges
    mutate(pathway_edge = melt(res$graph)$value) %>%
    dplyr::filter(pathway_edge == 1) %>%
    dplyr::select(-pathway_edge) %>%
    
    # compute additional properties
    purrr::pmap_dfr(function(
      source, target,
      dce, dce_stderr, dce_pvalue,
      cor, cor_pvalue, pcor, pcor_pvalue
    ) {
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
        cor = cor, cor_pvalue = cor_pvalue, pcor = pcor, pcor_pvalue = pcor_pvalue,
        distance = dist
      )
    }) %>%
    
    # format output
    arrange(dce_pvalue)
  
  # save results as dataframe
  df_final %>%
    write_csv(file.path(out.dir, glue::glue("dce_list_{appendix}.csv")))
  
  # network plot
  plot(
    res,
    labelsize = 0,
    node_color = "grey",
    node_border_size = 0.2,
    arrow_size = 0.02,
    highlighted_nodes = strsplit(perturbed.gene, ",")[[1]]
  ) +
    ggtitle(glue::glue("{pathway}: {perturbed.gene}"))
  ggsave(file.path(out.dir, glue::glue("network_{appendix}.pdf")), width = 6, height = 6)
  
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
      lab = df_volcano$edge, selectLab = dplyr::filter(df_volcano, df_volcano$affected)$edge,
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
    dplyr::filter(!is.infinite(distance)) %>%
    drop_na %>%
    ggplot(aes(x = -log10(dce_pvalue), y = distance, color = dce)) +
    geom_point() +
    xlab("-log10(DCE p-value)") +
    ylab("Graph distance perturbed-gene to DCE-edge") +
    ggtitle(glue::glue("{pathway}: {perturbed.gene}")) +
    theme_minimal()
  ggsave(file.path(out.dir, glue::glue("network_distances_{appendix}.pdf")), width = 8, height = 6)
  
}
