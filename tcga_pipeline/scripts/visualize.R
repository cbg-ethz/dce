library(tidyverse)


# load data
res <- readRDS(file=snakemake@input$dce_fname)
df.genes <- read_csv(snakemake@input$geneid_fname)

geneid.map <- setNames(as.character(df.genes$SYMBOL), df.genes$ENSEMBL)

# plot
p.list <- lapply(res, dce::plot.dce, nodename.map=geneid.map)

p <- cowplot::plot_grid(
  plotlist=p.list, labels=names(p.list),
  ncol=3
)
cowplot::save_plot(
  snakemake@output$plot_fname,
  p, ncol=3,
  base_height=10, base_asp=.8
)
