library(tidyverse)

library(KEGGgraph) # else: "data set `KEGGEdgeSubtype` not found"
library(Rgraphviz) # for `layoutGraph`

library(org.Hs.eg.db)


# load KEGG pathway
pw.id <- snakemake@wildcards$pathway
fname <- snakemake@output$xml_fname

if (!file.exists(fname)) {
  parts <- str_match(pw.id, "(\\D+)(\\d+)")
  KEGGgraph::retrieveKGML(pathwayid=pw.id, organism=parts[[2]], fname)
}


# set genesOnly=FALSE and do propagation and removal manually
# to keep genes which have only edges to non-genes
# (e.g. PTEN in BreastCancer pathway)
graph <- KEGGgraph::parseKGML2Graph(fname, genesOnly=FALSE, expandGenes=TRUE)
graph <- dce::propagate_gene_edges(graph)
graph


# plot pathway
p1 <- ggplotify::as.ggplot(~KEGGgraph::plotKEGGgraph(graph))
p2 <- ggplotify::as.ggplot(~KEGGgraph::KEGGgraphLegend())

p <- cowplot::plot_grid(p1, p2)
cowplot::save_plot(
  snakemake@output$plot_fname, p,
  ncol=2, nrow=1,
  base_height=8, base_width=10)


# relabel nodes (KEGG_ID -> Ensembl)
nodes.entrez <- KEGGgraph::translateKEGGID2GeneID(nodes(graph))

df.tmp <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys=nodes.entrez, keytype="ENTREZID",
  columns=c("ENTREZID", "ENSEMBL", "SYMBOL")
) %>%
  distinct(ENTREZID, .keep_all=TRUE)

stopifnot(length(nodes(graph)) == dim(df.tmp)[[1]])
df.tmp %>% write_csv(snakemake@output$geneid_fname)

df.tmp$ENSEMBL[is.na(df.tmp$ENSEMBL)] <- paste0("undef_entrez", df.tmp$ENTREZID[is.na(df.tmp$ENSEMBL)])
nodes.ensembl <- df.tmp %>% pull(ENSEMBL)

undef.count <- df.tmp %>% tally(grepl("undef", ENSEMBL)) %>% pull(n)
undef.count / dim(df.tmp)[1]

nodes(graph) <- nodes.ensembl


# save graph
as.data.frame(igraph::get.edgelist(igraph::igraph.from.graphNEL(graph))) %>%
  dplyr::rename(source=V1, sink=V2) %>%
  write_csv(snakemake@output$network_fname)
