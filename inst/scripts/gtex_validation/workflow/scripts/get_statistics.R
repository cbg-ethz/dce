###
# Compute performance measures (correlation, mean squared error).
###


print(snakemake@input[[1]])

mse_with_deconf = c()
cor_with_deconf = c()
mse_no_deconf = c()
cor_no_deconf = c()

for(file in snakemake@input){
  out = mget(load(file, envir=(NE. <- new.env())), envir=NE.)
  x1 = out$normal_fit_with_deconf$dce_pvalue
  x2 = out$extended_fit_with_deconf$dce_pvalue
  shared_genes = intersect(colnames(x1), colnames(x2))
  x1 = x1[shared_genes, shared_genes]
  x2 = x2[shared_genes, shared_genes]
  x1 = -log(c(x1)[is.na(c(x1))==FALSE] + 1e-10)
  x2 = -log(c(x2)[is.na(c(x2))==FALSE] + 1e-10)
  mse_with_deconf = c(mse_with_deconf, mean((x1-x2)^2))
  cor_with_deconf = c(cor_with_deconf, cor(x1, x2))

  y1 = out$normal_fit_no_deconf$dce_pvalue
  y2 = out$extended_fit_no_deconf$dce_pvalue
  shared_genes = intersect(colnames(y1), colnames(y2))
  y1 = y1[shared_genes, shared_genes]
  y2 = y2[shared_genes, shared_genes]
  y1 = -log(c(y1)[is.na(c(y1))==FALSE] + 1e-10)
  y2 = -log(c(y2)[is.na(c(y2))==FALSE] + 1e-10)
  mse_no_deconf = c(mse_no_deconf, mean((y1-y2)^2))
  cor_no_deconf = c(cor_no_deconf, cor(y1, y2))
}

save(mse_with_deconf,
     cor_with_deconf,
     mse_no_deconf,
     cor_no_deconf,
     file=snakemake@output[[1]])
