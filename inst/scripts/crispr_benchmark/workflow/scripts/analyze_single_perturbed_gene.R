library(tidyverse)


# gather parameters
fname.expr.wt <- snakemake@input$count_wt_file
fname.expr.mt <- snakemake@input$count_mt_file

cur.gene <- snakemake@wildcards$gene

out.dir <- snakemake@output$out_dir
dir.create(out.dir, recursive = TRUE)


# expression histograms
df.expr.wt <- read_csv(fname.expr.wt) %>%
  column_to_rownames_wrap("...1") %>% as.matrix
df.expr.mt <- read_csv(fname.expr.mt) %>%
  column_to_rownames_wrap("...1") %>% as.matrix


# aggregate isoforms, etc
rownames(df.expr.wt) <- sub("[.-].*", "", rownames(df.expr.wt))
df.expr.wt <- rowsum(df.expr.wt, rownames(df.expr.wt))

rownames(df.expr.mt) <- sub("[.-].*", "", rownames(df.expr.mt))
df.expr.mt <- rowsum(df.expr.mt, rownames(df.expr.mt))

stopifnot(sum(grepl('\\.', rownames(df.expr.wt))) == 0)
stopifnot(sum(grepl('\\.', rownames(df.expr.mt))) == 0)


# normalize gene counts
norm.method <- "raw" # c("raw", "tpm", "deseq2")

if (norm.method == "raw") {
  # nothing to do
} else if (norm.method == "tpm") {
  # retrieve gene lengths
  mart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
  df.lens <- biomaRt::getBM(
    attributes = c("hgnc_symbol", "start_position", "end_position"),
    filters = "hgnc_symbol",
    values = rownames(df.expr.mt),
    mart = mart) %>%
    distinct(hgnc_symbol, .keep_all = TRUE)
  df.lens$length <- df.lens$end_position - df.lens$start_position
  df.lens %>% head

  print(glue::glue("{length(setdiff(rownames(df.expr.mt), df.lens$hgnc_symbol))}/{length(rownames(df.expr.mt))} genes without entry in biomart"))

  genes <- intersect(df.lens$hgnc_symbol, rownames(df.expr.mt))
  gene.lengths <- unlist(split(df.lens$length, df.lens$hgnc_symbol))

  # compute TPMs
  compute.tpm <- function(counts, lens) {
    len.kb <- lens / 1000
    rpk <- counts / len.kb
    pm.scale <- colSums(rpk) / 1e6
    return(t(t(rpk) / pm.scale))
  }

  df.expr.wt <- compute.tpm(df.expr.wt[genes, ] + 1, gene.lengths[genes])
  df.expr.mt <- compute.tpm(df.expr.mt[genes, ] + 1, gene.lengths[genes])
} else if (norm.method == "deseq2") {
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = cbind(df.expr.wt, df.expr.mt),
    colData = data.frame(condition = c(rep("WT", dim(df.expr.wt)[[2]]), rep("MT", dim(df.expr.mt)[[2]]))),
    design = ~ condition
  )
  dds <- DESeq2::DESeq(dds)
  all.counts <- DESeq2::counts(dds, normalized = TRUE)

  tmp.wt <- all.counts[, colnames(all.counts)[grepl("Ctrl", colnames(all.counts))]]
  tmp.mt <- all.counts[, colnames(all.counts)[grepl("^((?!Ctrl).)*$", colnames(all.counts), perl = TRUE)]]

  stopifnot(all(dim(df.expr.wt) == dim(tmp.wt)))
  stopifnot(all(dim(df.expr.mt) == dim(tmp.mt)))

  df.expr.wt <- tmp.wt
  df.expr.mt <- tmp.mt
} else {
  stop(glue::glue("Invalid count normalization type: {norm.method}"))
}


# perturbed gene
for (gene in strsplit(cur.gene, ",")[[1]]) {
  expr.wt <- df.expr.wt[gene, ]
  expr.mt <- df.expr.mt[gene, ]

  data.frame(
    type = c(rep("WT", length(expr.wt)), rep("MT", length(expr.mt))),
    expression = c(expr.wt, expr.mt)
  ) %>%
    ggplot(aes(x = expression, fill = type)) +
      geom_histogram() +
      scale_y_log10() +
      ggtitle(glue::glue("Gene: {gene}")) +
      theme_minimal()

  ggsave(file.path(out.dir, glue::glue("expression_histogram_{gene}.pdf")))
}


# all genes
expr.wt.all <- df.expr.wt %>% as.numeric
expr.mt.all <- df.expr.mt %>% as.numeric

data.frame(
  type = c(rep("WT", length(expr.wt.all)), rep("MT", length(expr.mt.all))),
  expression = c(expr.wt.all, expr.mt.all)
) %>%
  ggplot(aes(x = expression, fill = type)) +
  geom_histogram() +
  scale_y_log10() +
  ggtitle(glue::glue("All genes")) +
  theme_minimal()

ggsave(file.path(out.dir, glue::glue("expression_histogram_all_genes.pdf")))


# save expression information
data.frame(
  perturbed.gene = cur.gene,
  gene.mean.WT = mean(expr.wt), # is only last gene
  gene.mean.MT = mean(expr.mt), # is only last gene
  all.mean.WT = mean(expr.wt.all),
  all.mean.MT = mean(expr.mt.all)
) %>%
  write_csv(file.path(out.dir, "expression_stats.csv"))
