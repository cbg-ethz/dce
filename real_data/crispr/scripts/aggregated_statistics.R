library(tidyverse)


# parameters
fname <- snakemake@input$fname

out_dir <- snakemake@output$out_dir
dir.create(out_dir, recursive = TRUE)


# read data
df_all <- read_csv(fname)
