###
# Compute differentially expressed genes using edgeR.
###


library(edgeR)
library(harmonicmeanp)

dge_net <- function(X,Y,G,stats=NULL) {
  if (is.null(stats)) {
    group <- factor(c(rep(1, nrow(X)), rep(2,nrow(Y))))
    y <- DGEList(counts=cbind(t(X), t(Y)), group=group)
    y <- calcNormFactors(y)
    design <- model.matrix( ~ group)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    stats <- glmQLFTest(fit, coef = 2)
  }
  fc <- stats$table$logFC
  pval <- stats$table$PValue
  edges <- which(G != 0, arr.ind = TRUE)
  dce <- dcep <- G
  for (i in seq_len(nrow(edges))) {
    x <- edges[i,1]
    y <- edges[i,2]
    dce[x,y] <- fc[y]-fc[x]
    dcep[x,y] <- min(pval[c(x,y)])
  }
  return(list(dce=-dce,dce_pvalue=dcep))
}
