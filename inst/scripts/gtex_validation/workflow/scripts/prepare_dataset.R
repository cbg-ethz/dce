###
# Preprocess GTEx datasets.
###


load(snakemake@input$gencode) #loads variable gencode

library(hash)
ens_to_symbol <- new.env(hash = TRUE)
for(i in 1:nrow(gencode)){
  ens_to_symbol[[gencode$gene_id[i]]] <- gencode$gene_name[i]
}

C <- read.csv(snakemake@input$covariates, header=T, sep="\t")
#if(dim(C)[2] < 300) {
#  next()
#}
all_data <- read.csv(snakemake@input$expressions, header = T, sep="\t")

expr <- t(all_data[,5:ncol(all_data)])
colnames(expr) = all_data$gene_id
for(i in 1:ncol(expr)){
  colnames(expr)[i] = ens_to_symbol[[colnames(expr)[i]]]
}
chrs <- all_data$X.chr

covariates <- t(C[,2:ncol(C)])

expr.normal <- expr
expr.unconfounded <- lm(expr ~ covariates)$residuals
chromosomes <- chrs

save(covariates, expr.normal, expr.unconfounded, chromosomes, file = snakemake@output[[1]])
