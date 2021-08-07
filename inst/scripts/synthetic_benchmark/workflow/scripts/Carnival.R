library(CARNIVAL)
library(edgeR)

carWrap <- function(X,Y,G,stats = NULL) {
  if (is.null(stats)) {
    group <- factor(c(rep(1, nrow(X)), rep(2,nrow(Y))))
    y <- DGEList(counts=cbind(t(X), t(Y)), group=group)
    y <- calcNormFactors(y)
    design <- model.matrix( ~ group)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    stats <- glmQLFTest(fit, coef = 2)
  }
  rownames(G) <- colnames(G) <- paste0("Node",1:ncol(G))
  Gmat2list <- function(G) {
    edges <- which(abs(G)>0,arr.ind=TRUE)
    L <- data.frame(source = rownames(G)[edges[,1]], interaction = rep(1,nrow(edges)), target = rownames(G)[edges[,2]])
    return(as_tibble(L))
  }
  carnivalOptions <- defaultLpSolveCarnivalOptions()
  L <- Gmat2list(G)
  data <- abs(stats$table$logFC)
  names(data) <- paste0("Node",1:length(data))
  data <- data[names(data) %in% unlist(c(L[,1],L[,3]))]
  res <- runInverseCarnival(measurements = data,
                            priorKnowledgeNetwork = L,
                            carnivalOptions = carnivalOptions)
  dce <- G*0
  idx <- which(res$weightedSIF$Node1 != "Perturbation")
  if (length(idx) > 1) {
    dce[cbind(res$weightedSIF$Node1,res$weightedSIF$Node2)[idx,]] <- res$weightedSIF$Weight[idx]*res$weightedSIF$Sign[idx]
  } else {
    idx1 <- which(res$weightedSIF$Node1 != "Perturbation")
    dce[res$weightedSIF$Node1[idx1],res$weightedSIF$Node2[idx1]] <- 100
  }
  pval <- 1-abs(dce)/100
  system("rm parsedData* lpFile*")
  return(list(dce=-dce,dce_pvalue=pval))
}
