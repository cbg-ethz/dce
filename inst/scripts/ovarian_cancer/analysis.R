###
# Download data, process sequencing data and compute DCEs.
###


library(tidyverse)

# download sequencing data
accession <- "GSE9891"
gse <- GEOquery::getGEO(accession)
files <- GEOquery::getGEOSuppFiles(accession)

system(glue::glue("cd {accession} && tar xf {accession}_RAW.tar"))

# extract meta information
sample.map <- Biobase::phenoData(gse[[1]])@data %>%
  rownames_to_column("sample") %>%
  dplyr::select(sample, title) %>%
  mutate(ID = as.numeric(str_remove(title, "^X"))) %>%
  dplyr::select(-title)

gene.map <- Biobase::featureData(gse[[1]])@data %>%
  as.data.frame %>%
  dplyr::select(ID, `Gene Symbol`) %>%
  mutate(symbol = str_remove(`Gene Symbol`, " *///.*"))

# process sequencing data
fns <- affy::list.celfiles(accession)
length(fns)
ab <- affy::ReadAffy(filenames = fns, celfile.path = accession)
ev <- vsn::vsnrma(ab)
data <- Biobase::exprs(ev)

# adjust row/column labels
colnames(data) <- str_remove(colnames(data), "\\.CEL\\.gz$")

tmp <- setNames(as.character(gene.map$symbol), gene.map$ID)
rownames(data) <- unname(tmp[rownames(data)])

# read covariates (from RetrieveCovariates.ipynb)
df.meta <- read_csv("covariates.csv") %>%
  inner_join(sample.map, "ID")

df.meta %>%
  group_by(Group) %>%
  tally

# stratify groups
group1 <- df.meta %>%
  dplyr::filter(Group == 1) %>%
  pull(sample)

group2 <- df.meta %>%
  dplyr::filter(Group %in% c(2,3,4,5,6)) %>%
  pull(sample)

# extract per-group counts
trans <- function(x) {
  floor(exp(x))
}

X.group1 <- data[, group1] %>%
  trans %>%
  t %>%
  as.data.frame

X.group2 <- data[, group2] %>%
  trans %>%
  t %>%
  as.data.frame

# read pathway
pathway <- igraph::graph_from_data_frame(read_csv("hsa04350.csv"))
nodes <- igraph::vertex.attributes(pathway)$name

# compute DCEs
common.genes <- intersect(nodes, rownames(data))

res <- dce::dce(
  igraph::induced_subgraph(pathway, common.genes),
  X.group1[, common.genes], X.group2[, common.genes],
  solver = "lm",
  lib_size = FALSE, latent = 0
)

plot(res, use_symlog = FALSE, nodesize = 30, labelsize = 2)
ggsave("output.pdf", width = 20, height = 20)
