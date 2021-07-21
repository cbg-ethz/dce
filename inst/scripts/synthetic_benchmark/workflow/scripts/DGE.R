library(edgeR)
library(harmonicmeanp)

dge_net <- function(X,Y,G,mode=1) {
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
  edges <- which(G == 1, arr.ind = TRUE)
  dce <- dcep <- G
  for (i in 1:ncol(edges)) {
    x <- edges[i,1]
    y <- edges[i,2]
    dce[x,y] <- fc[y]-fc[x]
    dcep[x,y] <- p.hmp(pval[c(x,y)])
  }
  return(list(dce=dce,dce_pvalue=dcep))
}

# p <- 15
# n <- 100
# G <- matrix(sample(c(0,1),p*p,replace=TRUE,prob=c(0.8,0.2)),p,p)
# G[lower.tri(G)] <- 0
# diag(G) <- 0
# rownames(G) <- colnames(G) <- 1:nrow(G)
#
# X <- matrix(rnbinom(n = n*p, size = 100, mu = 1000),100,p)
# Y <- matrix(rnbinom(n = n*p, size = 100, mu = 1000),100,p)
# diff <- sample(1:10,2)
# Y[,diff] <- rnbinom(n = 2*n, size = 100, mu = 100)
#
# colnames(G) <- rownames(G) <- colnames(Y) <- colnames(X) <- 1:ncol(X)
#
# test <- dge_net(X,Y,G)
#
# test2 <- dce(G,X,Y,solver = "lm")
