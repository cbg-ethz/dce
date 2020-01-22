library(tidyverse)


# load data
load(file=snakemake@input$dce_fname)

# plot
p.list <- lapply(res, dce::plot.dce)

p <- cowplot::plot_grid(
  plotlist=p.list, labels=names(p.list),
  ncol=3
)
cowplot::save_plot(
  snakemake@output$plot_fname,
  p, ncol=3,
  base_height=10, base_asp=.8
)
