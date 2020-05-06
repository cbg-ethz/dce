compute.mse <- function(y_pred, y_true) {
  return(mean((y_pred - y_true)^2))
}

compute.prec <- function(x, y, do = "prec") {
  tp <- sum(y < 0.05 & x != 0, na.rm = TRUE)
  fp <- sum(y < 0.05 & x == 0, na.rm = TRUE)
  if (do == "prec") {
    z <- tp/(tp+fp)
  } else if (do == "rec") {
    fn <- sum(y >= 0.05 & x != 0, na.rm = TRUE)
    z <- tp/(tp+fn)
  }
  return(z)
}

compute.prauc <- function(x, y, seq = 1000, NAweight = 1) {
  prec <- rec <- numeric(seq)
  auc <- 0
  alphas <- c(seq(0,1,length.out=seq-1), 2)
  for (i in seq_len(seq)) {
    alpha <- alphas[i]
    tp <- sum(y < alpha & x != 0, na.rm = TRUE)
    fp <- sum(y < alpha & x == 0, na.rm = TRUE)
    fn <- sum(y >= alpha & x != 0, na.rm = TRUE)
    prec[i] <- tp/(tp+fp)
    if (is.na(prec[i])) { prec[i] <- NAweight }
    rec[i] <- tp/(tp+fn)
    if (i > 1) {
      auc <- auc + (rec[i]-rec[i-1])*((prec[i]+prec[i-1])/2)
    }
  }
  return(auc)
}

compute.rocauc <- function(x, y, seq = 1000) {
  sens <- spec <- numeric(seq)
  auc <- 0
  alphas <- c(seq(0,1,length.out=seq-1), 2)
  for (i in seq_len(seq)) {
    alpha <- alphas[i]
    tp <- sum(y < alpha & x != 0, na.rm = TRUE)
    fp <- sum(y < alpha & x == 0, na.rm = TRUE)
    tn <- sum(y >= alpha & x == 0, na.rm = TRUE)
    fn <- sum(y >= alpha & x != 0, na.rm = TRUE)
    sens[i] <- tp/(tp+fn)
    spec[i] <- tn/(tn+fp)
    if (i > 1) {
      auc <- auc + (spec[i-1]-spec[i])*((sens[i]+sens[i-1])/2)
    }
  }
  return(auc)
}

apply.performance.measure <- function(df, func, label, ...) {
  return(
    data.frame(
      cor=func(df$truth, df$cor, ...),
      pcor=func(df$truth, df$pcor, ...),
      dce=func(df$truth, df$dce, ...),
      dce.lr=func(df$truth, df$dce.lr, ...),
      rand=func(df$truth, df$rand, ...)
    ) %>%
      mutate(type=label)
  )
}
