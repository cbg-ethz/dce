library(tidyverse)

library(TCGAbiolinks)
library(KEGGgraph)
library(igraph)
library(Rgraphviz) # for `layoutGraph`

library(cowplot)
library(ggplotify)

library(AnnotationDbi)
library(org.Hs.eg.db)

# load KEGG pathway
fname <- snakemake@output[["xml_file"]]
if (!file.exists(fname)) {
  parts <- str_match(fname, ".*/(\\D+)(\\d+).xml")
  KEGGgraph::retrieveKGML(pathwayid=parts[[3]], organism=parts[[2]], fname)
}

gR <- KEGGgraph::parseKGML2Graph(fname, genesOnly=TRUE, expandGenes=TRUE)

# plot pathway
p1 <- ggplotify::as.ggplot(~KEGGgraph::plotKEGGgraph(gR))
p2 <- ggplotify::as.ggplot(~KEGGgraph::KEGGgraphLegend())

p <- cowplot::plot_grid(p1, p2)
cowplot::save_plot(
  snakemake@output[["plot_file"]], p,
  ncol=2, nrow=1,
  base_height=8, base_width=10)

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
  write_csv(snakemake@output[["network_file"]])
