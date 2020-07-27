devtools::load_all("../../../dce")

library(tidyverse)

library(KEGGgraph)
library(org.Hs.eg.db)


# gather all human pathway IDs
pw.data <- KEGGREST::keggList("pathway", "hsa")
pw.ids <- sapply(names(pw.data), function(x) { strsplit(x, ":")[[1]][[2]] }, USE.NAMES = FALSE)

# TODO: do this in a more elegant way
pw.ids <- snakemake@params$pathways

# download, process and store all of them
target.dir <- strsplit(dirname(snakemake@output$graph_files[[1]]), "/")[[1]][[1]]
dir.create(file.path(target.dir, "kgml_files"), recursive = TRUE)

for (pw in pw.ids) {
  print(pw)

  # download
  fname.kgml <- file.path(target.dir, "kgml_files", glue::glue("{pw}.kgml"))
  if (!file.exists(fname.kgml)) {
    KEGGgraph::retrieveKGML(pathwayid = pw, organism = "hsa", fname.kgml)
  }

  # create network
  graph <- KEGGgraph::parseKGML2Graph(fname.kgml, genesOnly = FALSE, expandGenes = TRUE)
  graph <- dce::propagate_gene_edges(graph)

  # map gene IDs
  nodes.entrez <- KEGGgraph::translateKEGGID2GeneID(nodes(graph))

  df.tmp <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = nodes.entrez, keytype = "ENTREZID",
    columns = c("ENTREZID", "ENSEMBL", "SYMBOL")
  ) %>%
    distinct(ENTREZID, .keep_all=TRUE) %>%
    mutate(KEGG = nodes(graph))

  # identify nodes where KEGG ID has no SYMBOL mapping
  df.tmp$SYMBOL[is.na(df.tmp$SYMBOL)] <- paste0("undef_entrez", df.tmp$ENTREZID[is.na(df.tmp$SYMBOL)])

  undef.count <- df.tmp %>% tally(grepl("undef", SYMBOL)) %>% pull(n)
  print(glue::glue("Fraction of KEGG nodes (genes) without Gene Symbol: {undef.count / dim(df.tmp)[1]}"))

  # identify nodes where different KEGG IDs map to same Gene Symbol
  df.mapping.overlaps <- df.tmp %>%
    group_by(SYMBOL) %>%
    dplyr::filter(n() > 1)

  if (dim(df.mapping.overlaps)[[1]] > 0) {
    print("Multiple KEGG IDs map to same Gene Symbol:")
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
  geneid.map <- setNames(as.character(df.tmp$SYMBOL), df.tmp$KEGG)
  nodes(graph) <- unname(geneid.map[nodes(graph)])

  # save graph
  fname.csv <- file.path(target.dir, "csv_files", glue::glue("{pw}.csv"))
  as.data.frame(igraph::get.edgelist(igraph::igraph.from.graphNEL(graph))) %>%
    dplyr::rename(source=V1, sink=V2) %>%
    write_csv(fname.csv)
}
