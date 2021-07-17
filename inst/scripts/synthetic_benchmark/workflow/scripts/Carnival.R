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
  return(res)
}

G <- matrix(sample(c(0,1),10*10,replace=TRUE,prob=c(0.8,0.2)),10,10)
G[lower.tri(G)] <- 0
diag(G) <- 0
rownames(G) <- colnames(G) <- 1:nrow(G)

X <- matrix(rnbinom(n = 1000, size = 100, mu = 1000),100,10)
Y <- matrix(rnbinom(n = 1000, size = 100, mu = 1000),100,10)
diff <- sample(1:10,2)
Y[,diff] <- rnbinom(n = 200, size = 100, mu = 100)

test <- carWrap(X,Y,G)

