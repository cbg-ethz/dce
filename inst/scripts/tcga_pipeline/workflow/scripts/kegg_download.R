###
# Download pathways from KEGG.
###


library(tidyverse)

library(KEGGgraph) # else: "data set `KEGGEdgeSubtype` not found"
library(Rgraphviz) # for `layoutGraph`

library(org.Hs.eg.db)

devtools::load_all("../../../")


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


# plot pathway
p1 <- ggplotify::as.ggplot(~KEGGgraph::plotKEGGgraph(graph))
p2 <- ggplotify::as.ggplot(~KEGGgraph::KEGGgraphLegend())

p <- cowplot::plot_grid(p1, p2)
cowplot::save_plot(
  snakemake@output$plot_fname, p,
  ncol=2, nrow=1,
  base_height=8, base_width=10)


# gene id conversions (KEGG_ID -> Ensembl)
nodes.entrez <- KEGGgraph::translateKEGGID2GeneID(nodes(graph))

df.tmp <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys=nodes.entrez, keytype="ENTREZID",
  columns=c("ENTREZID", "ENSEMBL", "SYMBOL")
) %>%
  distinct(ENTREZID, .keep_all=TRUE) %>%
  mutate(KEGG=nodes(graph))

stopifnot(length(nodes(graph)) == dim(df.tmp)[[1]])
df.tmp %>% write_csv(snakemake@output$geneid_fname)


# identify nodes where KEGG ID has no ENSEMBL mapping
df.tmp$ENSEMBL[is.na(df.tmp$ENSEMBL)] <- paste0("undef_entrez", df.tmp$ENTREZID[is.na(df.tmp$ENSEMBL)])

undef.count <- df.tmp %>% tally(grepl("undef", ENSEMBL)) %>% pull(n)
print(glue::glue("Fraction of KEGG nodes (genes) without Ensembl ID: {undef.count / dim(df.tmp)[1]}"))


# identify nodes where different KEGG IDs map to same Ensembl ID
df.mapping.overlaps <- df.tmp %>%
  group_by(ENSEMBL) %>%
  dplyr::filter(n() > 1)

if (dim(df.mapping.overlaps)[[1]] > 0) {
  print("Multiple KEGG IDs map to same Ensembl ID:")
  print(df.mapping.overlaps)

  contraction.map <- df.mapping.overlaps %>%
    group_map(function(df.cur, name) {
      contract.idx <- match(df.cur$KEGG, nodes(graph))
      list(source.idx = contract.idx, target.idx = contract.idx[[1]])
    })
  contraction.mapping <- purrr::map(seq(1, length(nodes(graph))), function(x) {
    for (e in contraction.map) {
      if (x %in% e$source.idx) {
        return(e$target.idx)
      }
      return(x)
    }
  }) %>%
    unlist

  graph.ig <- igraph::igraph.from.graphNEL(graph)
  g.con <- igraph::contract(graph.ig, contraction.mapping, vertex.attr.comb = "first")

  # contraction leads to multiple edges between two nodes. Remove them
  g.con <- igraph::simplify(g.con)

  # contraction places isolated node with empty name in graph. Remove it
  isolated.nodes = which(igraph::degree(g.con) == 0)
  g.con <- igraph::delete.vertices(g.con, isolated.nodes)

  # contraction saves node names as list instead of vector. To retain names when converting to graphNEL, they have to be a vector
  raw.nodes <- unname(unlist(igraph::get.vertex.attribute(g.con, "name")))
  g.con <- igraph::set.vertex.attribute(igraph::delete_vertex_attr(g.con, "name"), "name", value = raw.nodes)

  graph <- igraph::igraph.to.graphNEL(g.con)
}


# relabel nodes
geneid.map <- setNames(as.character(df.tmp$ENSEMBL), df.tmp$KEGG)
nodes(graph) <- unname(geneid.map[nodes(graph)])


# save graph
as.data.frame(igraph::get.edgelist(igraph::igraph.from.graphNEL(graph))) %>%
  dplyr::rename(source=V1, sink=V2) %>%
  write_csv(snakemake@output$network_fname)
