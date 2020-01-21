library(tidyverse)


# load data
load(file=snakemake@input$dce_fname)

# plot
p <- dce::plot.dce(res)

cowplot::plot_grid(p, p) +
  ggsave(snakemake@output$plot_fname)
