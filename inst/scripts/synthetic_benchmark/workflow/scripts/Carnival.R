library(CARNIVAL)
library(edgeR)

carWrap <- function(X,Y,G) {
  Gmat2list <- function(G) {
    edges <- which(G==1,arr.ind=TRUE)
    L <- data.frame(source = rownames(G)[edges[,1]], interaction = rep(1,nrow(edges)), target = rownames(G)[edges[,2]])
  }
  carnivalOptions <- defaultLpSolveCarnivalOptions()
  group <- factor(c(rep(1,nrow(X)),rep(2,nrow(Y))))
  y <- DGEList(counts=cbind(t(X),t(Y)),group=group)
  y <- calcNormFactors(y)
  design <- model.matrix(~group)
  y <- estimateDisp(y,design)
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  L <- Gmat2list(G)
  data <- abs(qlf$table$logFC)
  names(data) <- rownames(G)
  res <- runInverseCarnival(measurements = data,
                            priorKnowledgeNetwork = L, 
                            carnivalOptions = carnivalOptions)
  dce <- G*0
  idx <- which(res$weightedSIF$Node1 != "Perturbation")
  dce[cbind(res$weightedSIF$Node1,res$weightedSIF$Node2)[idx,]] <- res$weightedSIF$Weight[idx]*res$weightedSIF$Sign[idx]
  pval <- 1-abs(dce)/100
  system("rm parsedData* lpFile*")
  return(list(dce=dce,dce_pvalue=pval))
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
# test <- carWrap(X,Y,G)
