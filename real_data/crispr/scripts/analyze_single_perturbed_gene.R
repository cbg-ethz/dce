library(tidyverse)


# gather parameters
fname.expr.wt <- snakemake@input$count_wt_file
fname.expr.mt <- snakemake@input$count_mt_file

cur.gene <- snakemake@wildcards$gene

out.dir <- snakemake@output$out_dir

# expression histograms
df.expr.wt <- read_csv(fname.expr.wt)
df.expr.mt <- read_csv(fname.expr.mt)

expr.wt <- df.expr.wt %>%
  dplyr::filter(X1 == cur.gene) %>%
  dplyr::select(-X1) %>%
  unlist %>%
  as.numeric
expr.mt <- df.expr.mt %>%
  dplyr::filter(X1 == cur.gene) %>%
  dplyr::select(-X1) %>%
  unlist %>%
  as.numeric

data.frame(
  type = c(rep("WT", length(expr.wt)), rep("MT", length(expr.mt))),
  expression = c(expr.wt, expr.mt)
) %>%
  ggplot(aes(x = expression, fill = type)) +
    geom_histogram() +
    scale_y_log10() +
    ggtitle(glue::glue("Gene: {cur.gene}")) +
    theme_minimal()

ggsave(file.path(out.dir, 'expression_histogram.pdf'))
