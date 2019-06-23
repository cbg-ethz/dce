library(tidyverse)

library(TCGAbiolinks)
library(KEGGgraph)
library(igraph)

library(AnnotationDbi)
library(org.Hs.eg.db)

# load KEGG pathway
fname <- "kegg_data/hsa05219.xml"
if (!file.exists(fname)) {
  KEGGgraph::retrieveKGML("05219", "hsa", fname)
}

gR <- KEGGgraph::parseKGML2Graph(fname, genesOnly=TRUE, expandGenes=TRUE)

# relabel nodes
convert.id <- function (x.kegg) {
  x.entrez <- KEGGgraph::translateKEGGID2GeneID(x.kegg)
  
  id.map <- AnnotationDbi::select(
    org.Hs.eg.db, 
    keys=x.entrez, keytype="ENTREZID",
    columns=c("ENTREZID", "ENSEMBL")
  ) %>%
    distinct(ENTREZID, .keep_all=TRUE) %>%
    column_to_rownames("ENTREZID")
  
  x.ensembl <- id.map[x.entrez,]
  return(x.ensembl)
}

nodes(gR) <- sapply(nodes(gR), convert.id)

# save graph
as.data.frame(igraph::get.edgelist(igraph::igraph.from.graphNEL(gR))) %>%
  dplyr::rename(source=V1, sink=V2) %>%
  write_csv("kegg_data/hsa05219.edgelist.csv")
