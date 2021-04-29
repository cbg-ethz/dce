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

# create overview table
dce_table <- purrr::map_dfr(res, function(x) {
  x %>%
    as.data.frame %>%
    mutate(
      source = geneid.map[as.character(source)],
      target = geneid.map[as.character(target)],
      dce_pvalue_str = format.pval(dce_pvalue)
    ) %>%
    arrange(desc(abs(dce))) #%>%
    # dplyr::filter(dce_pvalue <= 0.05 & abs(dce) > 1)
}, .id = "condition") %>%
  arrange(dce_pvalue) %>%
  mutate(
    dce_symlog = symlog(dce),
    edge_name = paste0(source, "->", target)
  )

dce_table %>%
  head

# volcano plot summary
color_map <- ifelse(
  dce_table$condition == "stage i",
  "red",
  ifelse(
    dce_table$condition == "stage ii",
    "green",
    ifelse(
      dce_table$condition == "stage iii",
      "blue",
      "grey"
    )
  )
)
names(color_map)[color_map == "red"] <- "stage i"
names(color_map)[color_map == "green"] <- "stage ii"
names(color_map)[color_map == "blue"] <- "stage iii"
table(color_map)

edge_selection <- dce_table %>%
  dplyr::filter(
    (
      grepl("DLL3", edge_name) |
        grepl("FGF", edge_name)
    ) & (
      abs(dce_symlog) > 1
    ) |
      dce_pvalue < 1e-50
  ) %>%
  pull(edge_name)

EnhancedVolcano::EnhancedVolcano(
  dce_table,
  lab = dce_table$edge_name, selectLab = edge_selection,
  x = "dce_symlog", y = "dce_pvalue",
  colCustom = color_map,
  pCutoff = .05, FCcutoff = 1,
  drawConnectors = TRUE,
  title = NULL, subtitle = NULL, caption = NULL,
  xlab = bquote(~symlog(DCE)), ylab = bquote(~-Log[10]~italic(pvalue)),
  legendLabels = c("NS", "DCE", "p-value", "p-value and DCE"),
  axisLabSize = 30, pointSize = 5, labSize = 8, legendLabSize = 30, legendIconSize = 10,
)
ggsave(file.path(outdir, glue::glue("volcanoplot.pdf")), width = 20, height = 20)
