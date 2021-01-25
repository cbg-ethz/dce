## paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6659630/
library(RcppParallel)
library(FastGGM)

Diff.omega <- function(omega1,omega2,n1,n2){
  W = matrix(0,ncol(omega1),ncol(omega1))
  for (i in 1:ncol(omega1)) {
    for (j in 1:ncol(omega1)){
      W[i,j] = (omega1[i,j]-omega2[i,j])/sqrt(((1/n1)*(omega1[i,i]*omega1[j,j]+(omega1[i,j])^2))+((1/n2)*(omega2[i,i]*omega2[j,j]+(omega2[i,j])^2)))
    }
  }
  return(W)
}

FDR_control<-function(diff,p,alpha){
  II = diag(p)
  w.upper = which(upper.tri(II))
  z=seq(0,floor(2*sqrt(log(p))+1),by=0.01)
  s=length(z)
  rr=1
  k=1
  while(rr>alpha & k<s){
    rr=(p^2-p)*(1-pnorm(z[k]))/max(sum(abs(diff[w.upper])>=z[k]),1)
    k=k+1
  }
  k=k-1
  selected.list=which(abs(diff[w.upper])>=z[k])
  return(selected.list)
}

######################standardized statistics for difference omega################################
FastGGM_Diff <-function(case,control,alpha=0.05){
  alpha = alpha
  n1 = nrow(case)
  n2 = nrow(control)
  p = ncol(case)
  II = diag(p)
  w.upper = which(upper.tri(II))
  w.mat = which(upper.tri(II), arr.ind = TRUE)
  genepair = data.frame(gene1 = colnames(case)[w.mat[, 1]],gene2 = colnames(case)[w.mat[, 2]])
  
  fastggm_case <- FastGGM_Parallel(case)
  fastggm_control <- FastGGM_Parallel(control)
  
  fastggm_case_control = Diff.omega(fastggm_case$precision,fastggm_control$precision,n1,n2)
  pvalue_fastggm_Diff = 2*pnorm(-abs(fastggm_case_control[w.upper]))
  
  results = data.frame(genepair = genepair,
                       case_precision = fastggm_case$precision[w.upper], 
                       case_precision_p = fastggm_case$p_precision[w.upper],
                       case_partialCor = fastggm_case$partialCor[w.upper], 
                       case_partialCor_p = fastggm_case$p_partialCor[w.upper],
                       control_precision = fastggm_control$precision[w.upper],
                       control_precision_p = fastggm_control$p_precision[w.upper],
                       control_partialCor = fastggm_control$partialCor[w.upper], 
                       control_partialCor_p = fastggm_control$p_partialCor[w.upper],
                       W = fastggm_case_control[w.upper],
                       pvalue_fastggm_Diff = pvalue_fastggm_Diff,
                       FDR_fastggm_Diff = p.adjust(pvalue_fastggm_Diff,"BH")
                       )
  dce <- fastggm_case_control
  dce_pvalue <- dce*0
  dce_pvalue[upper.tri(dce_pvalue)] <- pvalue_fastggm_Diff
  
  
  ## M = max(fastggm_case_control^2);
  ## # Calculate p-value
  ## Global_p_value = 1-exp(-1/sqrt(8*3.14159)*exp(-(M-4*log(p)+log(log(p)))/2))
  
  
  ## FDR_selection = FDR_control(fastggm_case_control,p,alpha)
  
  return(list(
      # M=M, Global_p_value=Global_p_value, results = results, W = fastggm_case_control[w.upper],FDR_results = results[FDR_selection,],
      dce=dce,dce_pvalue=dce_pvalue)) 
  
}
