library(dce)

if(snakemake@wildcards[['tissue1']] == snakemake@wildcards[['tissue2']]){
  quit()
}

ds1 = mget(load(snakemake@input[[1]], envir=(NE. <- new.env())), envir=NE.)
ds2 = mget(load(snakemake@input[[2]], envir=(NE. <- new.env())), envir=NE.)
X1 = ds1$expr.normal
H1 = ds1$covariates
X2 = ds2$expr.normal
H2 = ds2$covariates

pathways <- get_pathways(pathway_list = list(kegg = c(snakemake@config[['pathway']])))
pathway <- pathways[[1]]$graph
shared_genes <- intersect(nodes(pathway), intersect(colnames(X1), colnames(X2)))
subgraph = subGraph(shared_genes, pathway)
extended = subgraph

q = dim(H1)[2]
for(i in 1:q){
  extended = addNode(node=paste0('H', i), object=extended)
  for(node in shared_genes){
    extended = addEdge(from=paste0('H', i), to=node, graph=extended)
  }
}
colnames(H1) = paste0('H', seq_len(q))
colnames(H2) = paste0('H', seq_len(q))

glue::glue(
  "Covered nodes: {length(shared_genes)}/{length(nodes(pathway))}"
) %>% print
X1 = X1[, colnames(X1) %in% shared_genes]
X2 = X2[, colnames(X2) %in% shared_genes]

normal_fit_no_deconf <- dce::dce(subgraph, X1, X2, deconfounding=FALSE, test='vcovHC')
normal_fit_with_deconf <- dce::dce(subgraph, X1, X2, deconfounding=snakemake@config[['deconfounding']], test='vcovHC')

extended_fit_no_deconf <- dce::dce(extended, cbind(X1, H1), cbind(X2, H2), deconfounding=FALSE, test='vcovHC')
extended_fit_with_deconf <- dce::dce(extended, cbind(X1, H1), cbind(X2, H2), deconfounding=snakemake@config[['deconfounding']], test='vcovHC')

save(normal_fit_no_deconf,
     normal_fit_with_deconf,
     extended_fit_no_deconf,
     extended_fit_with_deconf,
     file=snakemake@output[[1]])