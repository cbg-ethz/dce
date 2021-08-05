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

custom_plot <- function(...) {
  plot(...) +
    theme(
      legend.title = element_text(size = rel(1.5)),
      legend.text = element_text(size = rel(1.2))
    )
}

p.list <- lapply(
  res, custom_plot,
  nodename_map = geneid.map,
  node_color = "grey", labelsize = 0,
  node_border_size = 0.2, arrow_size = 0.01,
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
      base_height = 4, base_asp = 2, limitsize = FALSE
    )
  })

# aggregate plot
cowplot::save_plot(
  file.path(outdir, "all_stages.pdf"),
  p, ncol = length(p.list),
  base_height = 30, base_asp = 1, limitsize = FALSE
)

# create overview table
dce_table <- purrr::map_dfr(res, function(x) {
  x %>%
    as.data.frame %>%
    mutate(
      source = geneid.map[as.character(source)],
      target = geneid.map[as.character(target)],
      dce_pvalue_str = format.pval(dce_pvalue)
    ) %>%
    arrange(desc(abs(dce)))
}, .id = "condition") %>%
  arrange(dce_pvalue) %>%
  mutate(
    dce_symlog = symlog(dce),
    edge_name = paste0(source, "→", target)
  ) %>%
  drop_na

dce_table %>%
  head

dce_table %>%
  dim

dce_table %>%
  dplyr::filter(dce_pvalue > 0.05 & abs(dce) < 1) %>%
  dim

# volcano plot summary
# custom labels because `selectLab` behaves weirdly for non-unique labels
volcano_table <- dce_table %>%
  mutate(
    edge_label = case_when(
      abs(dce_symlog) > 4.2 ~ edge_name,
      dce_pvalue < 1e-20 ~ edge_name,
      abs(dce_symlog) > 1 & -log10(dce_pvalue) > 12 ~ edge_name,
      dce_pvalue < 1e-14 ~ edge_name,
      dce_symlog < -4 ~ edge_name,
      TRUE ~ ""
    )
  )

# volcano_table %>%
#   arrange(dce_pvalue) %>%
#   view

# we also order them to have a nice legend in the plot, but only the first few
# rows, because otherwise the point overlaps produce misleading colors
old_dim <- dim(volcano_table)

slice_offset <- 50
volcano_table <- bind_rows(
  volcano_table %>%
    slice(1:slice_offset) %>%
    arrange(condition),
  volcano_table %>%
    slice(slice_offset + 1:n())
)

stopifnot(old_dim == dim(volcano_table))
volcano_table %>%
  head

# determine color scale
color_map <- case_when(
  volcano_table$condition == "stage i" ~ "red",
  volcano_table$condition == "stage ii" ~ "green",
  volcano_table$condition == "stage iii" ~ "blue",
  TRUE ~ "grey"
)
names(color_map)[color_map == "red"] <- "stage i"
names(color_map)[color_map == "green"] <- "stage ii"
names(color_map)[color_map == "blue"] <- "stage iii"
table(color_map)

# create actual plot
EnhancedVolcano::EnhancedVolcano(
  volcano_table,
  lab = volcano_table$edge_label, selectLab = volcano_table$edge_label,
  x = "dce_symlog", y = "dce_pvalue",
  colCustom = color_map,
  pCutoff = .05, FCcutoff = 1,
  drawConnectors = TRUE,
  title = NULL, subtitle = NULL, caption = NULL,
  xlab = bquote(~symlog(DCE)), ylab = bquote(~-log[10]~italic(pvalue)),
  axisLabSize = 30, pointSize = 4, labSize = 8, legendLabSize = 30, legendIconSize = 10
)
ggsave(file.path(outdir, glue::glue("volcanoplot.pdf")), width = 20, height = 20)
