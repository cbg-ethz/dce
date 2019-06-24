library(tidyverse)
library(TCGAbiolinks)

# general info
TCGAbiolinks::getGDCInfo()

df.projects <- TCGAbiolinks::getGDCprojects()
df.projects %>%
  filter(grepl("TCGA", project_id)) %>%
  head(1)

# prepare query
# TODO: download all project in single job
project.list <- c(str_match(snakemake@output[["data_dir"]], ".*/(.*)/.*")[[2]])

query <- TCGAbiolinks::GDCquery(
  project=project.list,
  data.category="Transcriptome Profiling",
  data.type="Gene Expression Quantification",
  experimental.strategy="RNA-Seq",
  workflow.type="HTSeq - Counts")

df.query <- query$results[[1]]
table(df.query$tissue.definition)

# download data
TCGAbiolinks::GDCdownload(query, directory=snakemake@output[["data_dir"]])

# convert to expression matrix
mat <- TCGAbiolinks::GDCprepare(
  query,
  directory=snakemake@output[["data_dir"]],
  summarizedExperiment=FALSE)

# save to file
df.query %>%
  dplyr::select(cases, tissue.definition) %>%
  write_csv(snakemake@output[["classification_file"]])

mat %>%
  dplyr::rename(gene=X1) %>%
  write_csv(snakemake@output[["expression_file"]])
