load(snakemake@input$gencode) #loads variable gencode

library(hash)
ens_to_symbol <- new.env(hash = TRUE)
for(i in 1:nrow(gencode)){
  ens_to_symbol[[gencode$gene_id[i]]] <- gencode$gene_name[i]
}

expr <- t(read.csv(snakemake@input$expressions, header=T, sep="\t"))
for(i in 1:ncol(expr)){
  colnames(expr)[i] = ens_to_symbol[[colnames(expr)[i]]]
}

C <- read.csv(snakemake@input$covariates, header=T, sep="\t")
covariates <- t(C[,2:ncol(C)])

common = intersect(rownames(expr), rownames(covariates))
covariates = covariates[common, ]
expr = expr[common,]

expr.normal <- expr
expr.unconfounded <- lm(expr ~ covariates)$residuals
save(covariates, expr.normal, expr.unconfounded, file = snakemake@output[[1]])