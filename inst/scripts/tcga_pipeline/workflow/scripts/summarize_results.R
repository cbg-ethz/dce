###
# Concatenate result csv files from all runs.
###


library(tidyverse)


# read data
dce_fname_list <- snakemake@input$dce_fname_list

purrr::map_dfr(dce_fname_list, function(fname) {
  res <- readRDS(file=fname)

  parts <- strsplit(basename(fname), "\\.")[[1]]
  dataset <- parts[[1]]
  pathway <- parts[[2]]

  purrr::map_dfr(res, function(x) {
    data.frame(pathway.pvalue = x$pathway.pvalue)
  }, .id = "condition") %>%
    mutate(dataset = dataset, pathway = pathway)
}) %>%
  write_csv(snakemake@output$enrichment_fname)
