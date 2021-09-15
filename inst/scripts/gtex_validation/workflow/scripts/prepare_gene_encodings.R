library(rtracklayer)
gencode <- rtracklayer::import(snakemake@input[[1]])
gencode = as.data.frame(gencode)
save(gencode, file = snakemake@output[[1]])