library(tidyverse)
library(ggraph)

devtools::load_all("../../../")


# parameters
outdir <- snakemake@output$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# load data
res <- readRDS(file=snakemake@input$dce_fname)
df.genes <- read_csv(snakemake@input$geneid_fname)

geneid.map <- setNames(as.character(df.genes$SYMBOL), df.genes$ENSEMBL)
geneid.map <- geneid.map[which(duplicated(names(geneid.map)) == FALSE)]

# prepare plots
dce.abs.max <- max(sapply(res, function(x) { max(abs(x$dce), na.rm = TRUE) }))
custom.limits <- c(-dce.abs.max, dce.abs.max)

p.list <- lapply(
  res, plot,
  nodename_map = geneid.map,
  node_color = "grey", labelsize = 1,
  edgescale_limits = custom.limits, use_symlog = TRUE
)

# labels.list <- purrr::map2(names(p.list), res, ~ glue::glue("{.x} (Enrichment: {.y$pathway_pvalue})")) %>% unlist
labels.list <- names(p.list) %>% unlist

p <- cowplot::plot_grid(
  plotlist = p.list, labels = labels.list,
  label_size = 40, ncol = length(p.list)
)

# individual plots
purrr::transpose(list(plot = p.list, label = labels.list)) %>%
  purrr::walk(function(x) {
    label_clean <- x$label %>% str_replace(" ", "_")
    print(label_clean)

    fname <- file.path(outdir, glue::glue("{label_clean}.pdf"))
    cowplot::save_plot(
      fname,
      x$plot,
      base_height = 10, base_asp = 1, limitsize = FALSE
    )
  })

# aggregate plot
cowplot::save_plot(
  file.path(outdir, "all_stages.pdf"),
  p, ncol = length(p.list),
  base_height = 30, base_asp = 1, limitsize = FALSE
)
