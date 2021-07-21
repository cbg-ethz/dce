library(edgeR)
library(harmonicmeanp)

dge_net <- function(X,Y,G) {
  group <- factor(c(rep(1,nrow(X)),rep(2,nrow(Y))))
  y <- DGEList(counts=cbind(t(X),t(Y)),group=group)
  y <- calcNormFactors(y)
  design <- model.matrix(~group)
  y <- estimateDisp(y,design)
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  fc <- qlf$table$logFC
  pval <- qlf$table$PValue
  pval[pval < 1e-50] <- 1e-50
  edges <- which(G != 0, arr.ind = TRUE)
  dce <- dcep <- G
  for (i in 1:nrow(edges)) {
    x <- edges[i,1]
    y <- edges[i,2]
    dce[x,y] <- fc[y]-fc[x]
    dcep[x,y] <- p.hmp(pval[c(x,y)])
  }
  return(list(dce=-dce,dce_pvalue=dcep))
}
