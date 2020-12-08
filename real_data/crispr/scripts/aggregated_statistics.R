library(tidyverse)


# parameters
fname <- snakemake@input$fname

out_dir <- snakemake@output$out_dir
dir.create(out_dir, recursive = TRUE)


# read data
df_all <- read_csv(fname)

df_all %>%
  head


# network distance plot
df_all %>%
  filter(!is.infinite(distance)) %>%
  drop_na() %>%
ggplot(aes(x = as.factor(distance), y = abs(dce))) +
  geom_boxplot() +
  scale_y_log10() +
  xlab("Graph distance: perturbed-gene to DCE-edge") +
  ylab("abs(DCE) (only >0.5)") +
  ylim(0.5, max(abs(df_all$dce))) +
  facet_wrap(~ study) +
  theme_minimal()
ggsave(file.path(out_dir, glue::glue("network_distances.pdf")), width = 8, height = 6)