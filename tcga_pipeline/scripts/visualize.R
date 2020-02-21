library(tidyverse)

library(dce)


# load data
res <- readRDS(file=snakemake@input$dce_fname)
df.genes <- read_csv(snakemake@input$geneid_fname)

geneid.map <- setNames(as.character(df.genes$SYMBOL), df.genes$ENSEMBL)

# plot
dce.abs.max <- max(sapply(res, function(x) { max(abs(x$dce)) }))
custom.limits <- dce::symlog(c(-dce.abs.max, dce.abs.max))

p.list <- lapply(res, plot, nodename.map=geneid.map, edgescale.limits=custom.limits)

p <- cowplot::plot_grid(
  plotlist=p.list, labels=names(p.list),
  label_size=40, ncol=length(p.list)
)

cowplot::save_plot(
  snakemake@output$plot_fname,
  p, ncol=length(p.list),
  base_height=20, base_asp=1, limitsize=FALSE
)
