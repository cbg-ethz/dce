# note regarding precision, recall, f1-score computation: https://github.com/dice-group/gerbil/wiki/Precision,-Recall-and-F1-measure

get.classification.counts <- function(df, col, alpha = 0.05) {
  # df.classification <- df %>%
  #   mutate_at(vars(-"truth"), ~ as.numeric(. <= alpha)) %>%
  #   mutate(truth = as.numeric(truth != 0))
  #
  # tmp.func <- function(y, y.hat) {
  #   c(
  #     tp = sum(y.hat == 1 & y == 1, na.rm = TRUE),
  #     fp = sum(y.hat == 1 & y == 0, na.rm = TRUE),
  #     tn = sum(y.hat == 0 & y == 0, na.rm = TRUE),
  #     fn = sum(y.hat == 0 & y == 1, na.rm = TRUE)
  #   )
  # }
  #
  # df.classification %>%
  #   map_dfc(~ tmp.func(df.classification$truth, .)) %>%
  #   select(-truth) %>%
  #   mutate(name = c("tp", "fp", "tn", "fn")) %>%
  #   column_to_rownames("name")
  # the above code is really cool but unfortunately has poor performance

  y <- as.numeric(df$truth != 0)
  y.hat <- as.numeric(df[[col]] <= alpha)

  list(
    tp = sum(y.hat == 1 & y == 1, na.rm = TRUE),
    fp = sum(y.hat == 1 & y == 0, na.rm = TRUE),
    tn = sum(y.hat == 0 & y == 0, na.rm = TRUE),
    fn = sum(y.hat == 0 & y == 1, na.rm = TRUE)
  )
}

compute.mse <- function(df, col) {
  return(mean((df[[col]] - df$truth)^2, na.rm = TRUE))
}

compute.precision <- function(df, col, alpha = .05) {
  c <- get.classification.counts(df, col, alpha = alpha)

  if (c$tp == 0) {
    if (c$fp == 0 && c$fn == 0 & any(!is.na(df[[col]]))) {
      # "there are no results and we (correctly) identify none"
      return(1)
    }
    # avoid dividing by 0
    return(0)
  }

  return(c$tp / (c$tp + c$fp))
}

compute.recall <- function(df, col, alpha = .05) {
  c <- get.classification.counts(df, col, alpha = alpha)

  if (c$tp == 0) {
    if (c$fp == 0 && c$fn == 0 & any(!is.na(df[[col]]))) {
      # "there are no results and we (correctly) identify none"
      return(1)
    }
    # avoid dividing by 0
    return(0)
  }

  return(c$tp / (c$tp + c$fn))
}

compute.f1score <- function(df, col, alpha = .05) {
  # handle special case
  c <- get.classification.counts(df, col, alpha = alpha)
  if (c$tp == 0) {
    if (c$fp == 0 && c$fn == 0 & any(!is.na(df[[col]]))) {
      # "there are no results and we (correctly) identify none"
      return(1)
    }
    # avoid dividing by 0
    return(0)
  }

  # actual computation
  prec <- compute.precision(df, col, alpha = alpha)
  reca <- compute.recall(df, col, alpha = alpha)

  return(2 * prec * reca / (prec + reca))
}

compute.prauc_old <- function(df, col, seq = 1000, NAweight = 1) {
  prec <- rec <- numeric(seq)
  auc <- 0
  alphas <- c(seq(0,1,length.out=seq-1), 2)
  for (i in seq_len(seq)) {
    alpha <- alphas[i]

    prec[i] <- compute.precision(df, col, alpha = alpha)
    if (is.na(prec[i])) { prec[i] <- NAweight }

    rec[i] <- compute.recall(df, col, alpha = alpha)

    if (i > 1) {
      auc <- auc + (rec[i]-rec[i-1])*((prec[i]+prec[i-1])/2)
    }
  }
  return(auc)
}

compute.prauc <- function(df, col, seq = 1000, NAweight = 1) {
  prec <- rec <- numeric(seq)
  auc <- 0
  alphas <- c(seq(0,1,length.out=seq-1), 2)
  for (i in seq_len(seq)) {
    alpha <- alphas[i]

    c <- get.classification.counts(df, col, alpha = alpha)

    prec[i] <- c$tp / (c$tp + c$fp)
    if (is.na(prec[i])) { prec[i] <- NAweight }
    rec[i] <- c$tp / (c$tp + c$fn)

    if (i > 1) {
      auc <- auc + (rec[i]-rec[i-1])*((prec[i]+prec[i-1])/2)
    }
  }
  return(auc)
}

compute.rocauc <- function(df, col, seq = 1000) {
  sens <- spec <- numeric(seq)
  auc <- 0
  alphas <- c(seq(0,1,length.out=seq-1), 2)
  for (i in seq_len(seq)) {
    alpha <- alphas[i]

    c <- get.classification.counts(df, col, alpha = alpha)

    sens[i] <- c$tp / (c$tp + c$fn)
    spec[i] <- c$tn / (c$tn + c$fp)

    if (i > 1) {
      auc <- auc + (spec[i-1]-spec[i])*((sens[i]+sens[i-1])/2)
    }
  }
  return(auc)
}

apply.performance.measure <- function(df, methods, func, label, ...) {
  purrr::map_dfc(methods, function(method) {
    tmp <- data.frame(func(df, method, ...))
    names(tmp) <- method
    tmp
  }) %>%
    mutate(type=label)
}
