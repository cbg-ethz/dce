library(tidyverse)

library(dce)

library(gridExtra)

library(epiNEM)

# load data
res <- readRDS(file=snakemake@input$dce_fname)

df.genes <- read_csv(snakemake@input$geneid_fname)

geneid.map <- setNames(as.character(df.genes$SYMBOL), df.genes$ENSEMBL)

# plot

type<- "graph"
if (type == "graph") {
    dce.abs.max <- max(sapply(res, function(x) { max(abs(x$dce)) }))
    custom.limits <- dce::symlog(c(-dce.abs.max, dce.abs.max))
    
    p.list <- lapply(res, plot, nodename.map=geneid.map, edgescale.limits=custom.limits)
    
    p <- cowplot::plot_grid(
                      plotlist=p.list, labels=names(p.list),
                      label_size=40, ncol=length(p.list)
                  )
    
    cowplot::save_plot(
                 snakemake@output$plot_fname,
                 p, ncol=length(p.list),
                 base_height=20, base_asp=1, limitsize=FALSE
             )
} else if (type == "new") {
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

p <- list()
for (i in seq_len(length(res))) {
    x <- res[[i]]
    if (i == length(res)) {
        colorkey <- list(space="right")
    } else {
        colorkey <- NULL
    }
    p[[i]] <- dce::plotDce(x, type = "heatmap", log = TRUE, col = "RdBu",
                           bordercol = "transparent", aspect = "iso",
                           colorkey = colorkey, main = names(res)[i],
                           cexMain = 4, clusterx = res[[1]]$dce,
                           genelabels = geneid.map, scalefac = dcemax, par.settings=list(layout.heights=list(top.padding=-20)))
}
pdf(snakemake@output$heatmap_fname, width = 60, height = 20)
gridExtra::grid.arrange(grobs=p, ncol=length(res))
dev.off()
