library(tidyverse)
library(TCGAbiolinks)

# general info
TCGAbiolinks::getGDCInfo()

df.projects <- TCGAbiolinks::getGDCprojects()
df.projects %>%
  filter(grepl("TCGA", project_id)) %>%
  head(1)

# prepare query
project.list <- c("TCGA-BLCA")

query <- TCGAbiolinks::GDCquery(
  project=project.list,
  data.category="Transcriptome Profiling",
  data.type="Gene Expression Quantification",
  experimental.strategy="RNA-Seq",
  workflow.type="HTSeq - Counts")

df.query <- query$results[[1]]
table(df.query$tissue.definition)

# download data
download_path <- 'tcga_data'
TCGAbiolinks::GDCdownload(query, directory=paste(c(download_path, "GDCdata"), collapse="/"))

# convert to expression matrix
mat <- TCGAbiolinks::GDCprepare(
  query, 
  directory=paste(c(download_path, "GDCdata"), collapse="/"), 
  summarizedExperiment=FALSE)

# save to file
df.query %>%
  dplyr::select(cases, tissue.definition) %>%
  write_csv(paste(c(download_path, "case_classifications.csv"), collapse="/"))

mat %>%
  dplyr::rename(gene=X1) %>%
  write_csv(paste(c(download_path, "expression_matrix.csv"), collapse="/"))