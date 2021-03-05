library(tidyverse)
library(ggraph)

devtools::load_all("../../../")


# load data
res <- readRDS(file=snakemake@input$dce_fname)
df.genes <- read_csv(snakemake@input$geneid_fname)

geneid.map <- setNames(as.character(df.genes$SYMBOL), df.genes$ENSEMBL)
geneid.map <- geneid.map[which(duplicated(names(geneid.map)) == FALSE)]

# plot
dce.abs.max <- max(sapply(res, function(x) { max(abs(x$dce), na.rm = TRUE) }))
custom.limits <- c(-dce.abs.max, dce.abs.max)

p.list <- lapply(res, plot, nodename_map = geneid.map, edgescale_limits = custom.limits, use_symlog = FALSE)
labels.list <- purrr::map2(names(p.list), res, ~ glue::glue("{.x} (Enrichment: {.y$pathway.pvalue})")) %>% unlist

p <- cowplot::plot_grid(
  plotlist = p.list, labels = labels.list,
  label_size = 40, ncol = length(p.list)
)

cowplot::save_plot(
  snakemake@output$plot_fname,
  p, ncol = length(p.list),
  base_height = 30, base_asp = 1, limitsize = FALSE
)
