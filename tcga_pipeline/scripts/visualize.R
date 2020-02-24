library(tidyverse)

library(dce)


# load data
res <- readRDS(file=snakemake@input$dce_fname)
df.genes <- read_csv(snakemake@input$geneid_fname)

geneid.map <- setNames(as.character(df.genes$SYMBOL), df.genes$ENSEMBL)
geneid.map <- geneid.map[which(duplicated(names(geneid.map)) == FALSE)]

# plot
dce.abs.max <- max(sapply(res, function(x) { max(abs(x$dce)) }))
custom.limits <- c(-dce.abs.max, dce.abs.max)

p.list <- lapply(res, plot, nodename.map=geneid.map, edgescale.limits=custom.limits)
labels.list <- purrr::map2(names(p.list), res, ~ glue::glue("{.x} (Enrichment: {.y$pathway.pvalue})")) %>% unlist

p <- cowplot::plot_grid(
  plotlist=p.list, labels=labels.list,
  label_size=40, ncol=length(p.list)
)

cowplot::save_plot(
  snakemake@output$plot_fname,
  p, ncol=length(p.list),
  base_height=30, base_asp=1, limitsize=FALSE
)

# heatmap:
library(gridExtra)
library(epiNEM)

dcemax <- lapply(res, function(x) {
    tmp <- x$dce
    tmp <- as.vector(tmp[which(tmp != 0)])
    geq0 <- which(tmp > 0)
    leq0 <- which(tmp < 0)
    tmp[geq0] <- log(tmp[geq0] + 1)
    tmp[leq0] <- -log(-tmp[leq0] + 1)
    return(log = abs(tmp))
})
dcemax <- max(unlist(dcemax))

n <- length(res)

dcesum <- 0
for (i in seq_len(n)) {
    dcesum <- res[[i]]$dce + dcesum
}
dcesum[is.na(dcesum)] <- 0

pdf(snakemake@output$heatmap_fname, width = n*20, height = 7*n)
p <- list()
for (i in seq_len(n)) {
    x <- res[[i]]
    p[[i]] <- dce::plotDce(x, type = "heatmap", log = TRUE, col = "RdBu",
                           bordercol = "transparent", aspect = "iso",
                           colorkey = NULL, main = names(res)[i],
                           cexMain = n, clusterx = dcesum,
                           genelabels = geneid.map, scalefac = dcemax)
}
gridExtra::grid.arrange(grobs=p, ncol=length(res), widths = rep(22, 3))
dev.off()

## alternative:
if (FALSE) {
    genesymbol <- as.character(df.genes$SYMBOL)
    names(genesymbol) <- df.genes$ENSEMBL

    dcemax <- lapply(res, function(x) {
        tmp <- x$dce
        tmp <- as.vector(tmp[which(tmp != 0)])
        geq0 <- which(tmp > 0)
        leq0 <- which(tmp < 0)
        tmp[geq0] <- log(tmp[geq0] + 1)
        tmp[leq0] <- -log(-tmp[leq0] + 1)
        return(log = abs(tmp))
    })
    dcemax <- max(unlist(dcemax))
    
    pdf(snakemake@output$plot_fname, width = length(res)*50, height = 50)
    top <- NULL
    par(mfrow=c(1,length(res)))
    for (i in seq_len(length(res))) {
        dce::plotDce(res[[i]], log = TRUE, logfun = log10, scalefac = dcemax,
                     nodelabel = as.list(genesymbol), main = names(res)[i])
    }
    dev.off()
}

