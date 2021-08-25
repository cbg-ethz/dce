rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(rtracklayer)
gencode <- rtracklayer::import("./gencode.v26.GRCh38.genes.gtf")
gencode = as.data.frame(gencode)
library(hash)
ens_to_symbol <- new.env(hash = TRUE)
for(i in 1:nrow(gencode)){
  ens_to_symbol[[gencode$gene_id[i]]] <- gencode$gene_name[i]
}

fns <- list.files("./GTEx_Analysis_v8_sQTL_covariates/", full.names = T)
datasets <- list()
for(fn in fns) {
  C <- read.csv(fn, header=T, sep="\t")
  if(dim(C)[2] < 300) {
    next()
  }
  tissue_name <- strsplit(strsplit(fn, ".", fixed = T)[[1]][2], "/")[[1]][4]
  d_fn <- paste("./GTEx_Analysis_v8_eQTL_expression_matrices/", tissue_name, ".v8.normalized_expression.bed.gz", sep="")
  all_data <- read.csv(d_fn, header = T, sep="\t")
  
  expr <- t(all_data[,5:ncol(all_data)])
  colnames(expr) = all_data$gene_id
  for(i in 1:ncol(expr)){
    colnames(expr)[i] = ens_to_symbol[[colnames(expr)[i]]]
  }
  chrs <- all_data$X.chr
  
  covariates <- t(C[,2:ncol(C)])
  
  L <- list()
  L$covariates <- covariates
  L$expr.normal <- expr
  L$expr.unconfounded <- lm(expr ~ covariates)$residuals
  L$chromosomes <- chrs
  datasets[[tissue_name]] <- L
  print(tissue_name)
  print(dim(expr))
}

save(file = "all_tissues.preprocessed.data2.RData.gz", compress = T, gencode, datasets)
