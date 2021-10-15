###
# Download gene expression data from TCGA
# (https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga).
###


library(tidyverse)


# general info
TCGAbiolinks::getGDCInfo()

df.projects <- TCGAbiolinks::getGDCprojects()
df.projects %>%
  filter(grepl("TCGA", project_id)) %>%
  head(1)

# prepare query
project.list <- snakemake@config$projects

query <- TCGAbiolinks::GDCquery(
  project=project.list,
  data.category="Transcriptome Profiling",
  data.type="Gene Expression Quantification",
  experimental.strategy="RNA-Seq",
  workflow.type="HTSeq - Counts")

df.query <- query$results[[1]]
table(df.query$tissue.definition)

# download data
TCGAbiolinks::GDCdownload(query, directory=snakemake@output$data_dir)

# access data
data <- TCGAbiolinks::GDCprepare(query, directory=snakemake@output$data_dir)

# save to file
HDF5Array::saveHDF5SummarizedExperiment(data, dir=snakemake@output$hdf5_dname, replace=TRUE)
