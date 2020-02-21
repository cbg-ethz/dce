library(tidyverse)

library(dce)

library(gridExtra)

library(epiNEM)

# load data
res <- readRDS(file=snakemake@input$dce_fname)

df.genes <- read_csv(snakemake@input$geneid_fname)

geneid.map <- setNames(as.character(df.genes$SYMBOL), df.genes$ENSEMBL)

# plot
dce.abs.max <- max(sapply(res, function(x) { max(abs(x$dce)) }))
custom.limits <- c(-dce.abs.max, dce.abs.max)

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

## alternative:

## graph <- as(res[[1]]$graph, "matrix")

## genesymbol <- as.character(df.genes$SYMBOL)
## names(genesymbol) <- df.genes$ENSEMBL
## genesymbol <- genesymbol[which(names(genesymbol) %in% rownames(graph))]

## tmp <- NULL
## for (i in seq_len(length(res))) {
##     tmp <- cbind(tmp, NA, res[[i]]$dce)
## }
## tmp <- tmp[, -1]
## rownames(tmp) <- genesymbol
## colnames(tmp) <- c(genesymbol, "", genesymbol, "", genesymbol)
## colcode <- c(rep(rgb(0.5,0,0), nrow(tmp)), 1,
##              rep(rgb(0.75,0,0), nrow(tmp)), 1,
##              rep(rgb(1,0,0), nrow(tmp)))

## geq0 <- which(tmp > 0)
## leq0 <- which(tmp < 0)
## tmp[geq0] <- log(tmp[geq0] + 1)
## tmp[leq0] <- -log(-tmp[leq0] + 1)

## one heatmap:

## pdf(snakemake@output$plot_fname, width = 60, height = 20)
## epiNEM::HeatmapOP(tmp, Colv = FALSE, Rowv = FALSE, col = "RdBu", bordercol = "transparent", colSideColors = colcode, aspect = "iso")
## dev.off()

## 3 heatmaps:

## pdf(snakemake@output$plot_fname, width = 60, height = 20)
## p1 <- epiNEM::HeatmapOP(tmp[, 1:ncol(graph)], Colv = FALSE, Rowv = FALSE, col = "RdBu", bordercol = "transparent", aspect = "iso", colorkey = NULL)
## p2 <- epiNEM::HeatmapOP(tmp[, (ncol(graph)+2):(ncol(graph)*2+1)], Colv = FALSE, Rowv = FALSE, col = "RdBu", bordercol = "transparent", aspect = "iso", colorkey = NULL)
## p3 <- epiNEM::HeatmapOP(tmp[, (ncol(graph)*2+4):(ncol(graph)*3+2)], Colv = FALSE, Rowv = FALSE, col = "RdBu", bordercol = "transparent", aspect = "iso")
## gridExtra::grid.arrange(p1, p2, p3, ncol=3)
## dev.off()

## alternative:

## genesymbol <- as.character(df.genes$SYMBOL)
## names(genesymbol) <- df.genes$ENSEMBL
## genesymbol <- as.list(genesymbol)

## graph <- as(res[[1]]$graph, "matrix")
## graph1 <- which(graph == 1, arr.ind = TRUE)
## graph <- apply(graph1, 1, function(x) {
##     y <- paste0(rownames(graph)[x[1]], "=", colnames(graph)[x[2]])
##     return(y)
## })

## bordercol <- lapply(genesymbol, function(x) return("black"))
## nodecol <- lapply(genesymbol, function(x) return("transparent"))
## genesymbolFull <- lapply(genesymbol, function(x) return(""))

## dcemax <- lapply(res, function(x) {
##     tmp <- x$dce
##     tmp <- as.vector(tmp[which(tmp != 0)])
##     geq0 <- which(tmp > 0)
##     leq0 <- which(tmp < 0)
##     tmp[geq0] <- log(tmp[geq0] + 1)
##     tmp[leq0] <- -log(-tmp[leq0] + 1)
##     return(abs(tmp))
## })
## dcemax <- max(unlist(dcemax))

## pdf(snakemake@output$plot_fname, width = length(res)*50, height = 50)
## top <- NULL # seq_len(length(graph))
## par(mfrow=c(1,length(res)))
## for (i in seq_len(length(res))) {
##     tmp <- res[[i]]$dce
##     tmp <- as.vector(tmp[which(tmp != 0)])
##     geq0 <- which(tmp > 0)
##     leq0 <- which(tmp < 0)
##     tmp[geq0] <- log(tmp[geq0] + 1)
##     tmp[leq0] <- -log(-tmp[leq0] + 1)
##     tmpBlue <- tmp/dcemax
##     tmpRed <- -tmp/dcemax
##     tmpBlue[which(tmpBlue < 0)] <- 0
##     tmpRed[which(tmpRed < 0)] <- 0
##     edgecol <- rgb(tmpRed,0,tmpBlue,apply(cbind(tmpBlue+tmpRed, 0.1), 1, max))
##     edgewidth <- apply(cbind(tmpBlue, tmpRed), 1, max)
##     mnem::plotDnf(graph, edgecol = edgecol,
##                   edgewidth = (edgewidth*4)+1,
##                   bordercol = bordercol, nodecol = nodecol,
##                   nodelabel = genesymbol)
##     top <- unique(c(top, which(edgewidth > quantile(edgewidth, 0.99))))
## }
## dev.off()

## pdf(paste0(gsub(".pdf", "", snakemake@output$plot_fname), "_zoom.pdf"), width = length(res)*5, height = 5)
## par(mfrow=c(1,length(res)))
## for (i in seq_len(length(res))) {
##     tmp <- res[[i]]$dce
##     tmp <- as.vector(tmp[which(tmp != 0)])
##     tmpBlue <- tmp/dcemax
##     tmpRed <- -tmp/dcemax
##     tmpBlue[which(tmpBlue < 0)] <- 0
##     tmpRed[which(tmpRed < 0)] <- 0
##     edgecol <- rgb(tmpRed,0,tmpBlue,apply(cbind(tmpBlue+tmpRed, 0.1), 1, max))
##     edgewidth <- apply(cbind(tmpBlue, tmpRed), 1, max)
##     graphTop <- graph[top]
##     edgecol <- edgecol[top]
##     edgewidth <- edgewidth[top]
##     bordercol <- bordercol[top]
##     nodecol <- nodecol[top]
##     mnem::plotDnf(graphTop, edgecol = edgecol,
##                   edgewidth = (edgewidth*4)+1,
##                   bordercol = bordercol, nodecol = nodecol,
##                   nodelabel = genesymbol, edgelabel = round(tmp[top], -2), labelcol = rep("black", length(top)))
## }
## dev.off()
