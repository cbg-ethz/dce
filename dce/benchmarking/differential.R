library(tidyverse)
library(magrittr)
library(graph)

devtools::load_all("..")

# set.seed(42) # do not use this, if you append independent runs!


# parse commandline arguments
"
Benchmark DCE performance and runtime.

Usage:
  differential.R
  differential.R --variable NAME --values VALUES --append TRUE/FALSE --replicates NUMBER

Options:
  -h --help        Show this screen.
  --variable NAME  Which property to vary [default: node.num].
  --values VALUES  What values to assign to varied property [default: 20,50,100].
  --append FALSE   If TRUE appends the results of this/these run(s) to an existing results file.
  --replicates 100 Number of simulation runs.
" -> doc

arguments <- docopt::docopt(doc)


# global parameters
node.num <- 10
wt.samples <- 200
mt.samples <- 200
beta.magnitude <- 1
dispersion <- 1
adjustment.type <- "parents"
dist.mean <- 1000
sample.kegg <- FALSE
append <- FALSE
replicate.count <- 100
perturb <- 0
true.positives <- 0.5


# parse parameters
varied.parameter <- arguments$variable
parameter.list <- unlist(
  purrr::map(strsplit(arguments$values, ",")[[1]], type.convert)
)
replicate.count <- arguments$replicates
append <- as.logical(arguments$append)

print(glue::glue("Benchmark parameters:"))
print(glue::glue("  Varied parameter: {varied.parameter}"))
print(glue::glue("  Parameter: {parameter.list}"))


# helper functions
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

# do benchmarking
if (sample.kegg) {
    kegg.dag <- readRDS("pathways.rds")
    node.num <- 10^9
    replicate.count <- length(kegg.dag)
}

seed.list <- sample(seq_len(10^9), replicate.count)

df.bench <- purrr::pmap_dfr(
  list(parameter=rep(parameter.list, each=replicate.count), index=rep(seq_len(replicate.count), length(parameter.list))),
  purrr::possibly(
    function(parameter, index) {
      # handle parametrization
      switch(
        varied.parameter,
        node.num={ node.num <- parameter },
        wt.samples={ wt.samples <- parameter },
        mt.samples={ mt.samples <- parameter },
        beta.magnitude={ beta.magnitude <- parameter },
        dispersion={ dispersion <- parameter },
        adjustment.type={ adjustment.type <- parameter },
        perturb={ perturb <- parameter },
        true.positives={ true.positives <- parameter }
      )
      print(glue::glue("node.num={node.num} wt.samples={wt.samples} mt.samples={mt.samples} beta.magnitude={beta.magnitude} dispersion={dispersion} adjustment.type={adjustment.type} perturb={perturb}"))
      set.seed(seed.list[index])

      # create graphs
      edge.prob <- runif(1, 0, 1)

      negweight.range <- c(-beta.magnitude, 0)
      posweight.range <- c(0, beta.magnitude)

      if (sample.kegg) {
          wt.graph <- kegg.dag[[sample(1:length(kegg.dag), 1)]]
          kegg.order <- order(apply(wt.graph, 1, sum) - apply(wt.graph, 2, sum), decreasing = TRUE)
          wt.graph <- wt.graph[kegg.order, kegg.order]
          wt.graph[lower.tri(wt.graph)] <- 0
          diag(wt.graph) <- 0
          wt.graph <- wt.graph[seq_len(min(node.num, nrow(wt.graph))), seq_len(min(node.num, nrow(wt.graph)))]
          wt.graph <- as(wt.graph, "graphNEL")
          wt.graph <- resample_edge_weights(wt.graph, negweight.range, posweight.range)
      } else {
          wt.graph <- create_random_DAG(node.num, edge.prob, negweight.range, posweight.range)
          while(length(wt.graph@edgeData@data) <= 1) {
              wt.graph <- create_random_DAG(node.num, edge.prob, negweight.range, posweight.range)
          }
      }
      mt.graph <- resample_edge_weights(wt.graph, negweight.range, posweight.range, tp = true.positives)

      # generate data
      wt.X <- simulate_data(wt.graph, n=wt.samples, dist.dispersion=dispersion, dist.mean=dist.mean)
      mt.X <- simulate_data(mt.graph, n=mt.samples, dist.dispersion=dispersion, dist.mean=dist.mean)

      dispersion.estimate <- estimateTheta(rbind(wt.X, mt.X))
      mean.estimate <- mean(rbind(wt.X, mt.X))  

      # sanity checks
      if (any(is.nan(wt.X)) || any(is.nan(mt.X))) {
        stop("Malformed simulated data")
      }

      # perturb dag
      if (perturb > 0) {
          p.dag <- as(wt.graph, "matrix")
          p.dag[which(p.dag != 0)] <- 1
          candidates <- intersect(which(p.dag == 0),
                                  which(upper.tri(p.dag) == TRUE))
          turn <- sample(candidates, floor(perturb*length(candidates)))
          p.dag[turn] <- 1
          p.dag<- as(p.dag, "graphNEL")
      } else if (perturb < 0) {
          p.dag <- as(wt.graph, "matrix")
          p.dag[which(p.dag != 0)] <- 1
          candidates <- which(p.dag != 0)
          sample.n <- floor(abs(perturb)*length(candidates))
          if (length(candidates) - sample.n <= 1) {
              sample.n <- length(candidates)
          }
          turn <- sample(candidates, sample.n)
          p.dag[turn] <- 0
          p.dag <- as(p.dag, "graphNEL")
      } else {
          p.dag <- wt.graph
      }

      # compute graph features
      tmp <- as(wt.graph, "matrix")
      tmp[which(tmp != 0)] <- 1
      graph.density <- sum(tmp) / ((dim(tmp)[1] * (dim(tmp)[1] - 1)) / 2)

      # run models
      ground.truth <- list(dce=trueEffects(mt.graph) - trueEffects(wt.graph))

      time.tmp <- Sys.time()
      res.cor <- list(dce=cor(mt.X) - cor(wt.X))
      res.cor$dce.pvalue <- pcor_perm(wt.X, mt.X, fun = cor)  
      time.cor <- as.integer(difftime(Sys.time(), time.tmp, units="secs"))

      time.tmp <- Sys.time()
      res.pcor <- list(dce=pcor(mt.X) - pcor(wt.X))
      res.pcor$dce.pvalue <- pcor_perm(wt.X, mt.X, fun = pcor)  
      time.pcor <- as.integer(difftime(Sys.time(), time.tmp, units="secs"))

      time.tmp <- Sys.time()
      res.dce <- try(dce::dce.nb(
        p.dag, wt.X, mt.X,
        adjustment.type = adjustment.type
      ))
      if (length(grep("Error", res.dce)) > 0) {
          print("this seed produced an error:")
          print(seed.list[index])
          print("error msg:")
          print(res.dce)
      }
      time.dce <- as.integer(difftime(Sys.time(), time.tmp, units="secs"))

      time.tmp <- Sys.time()
      tmp <- as.matrix(ground.truth$dce)
      tmp[which(as.matrix(ground.truth$dce) != 0)] = (
        runif(sum(as.matrix(ground.truth$dce) != 0), negweight.range[1], posweight.range[2]) -
        runif(sum(as.matrix(ground.truth$dce) != 0), negweight.range[1], posweight.range[2])
      )
      tmp2 <- tmp*0
      tmp2[which(as.matrix(ground.truth$dce) != 0)] <- runif(sum(as.matrix(ground.truth$dce) != 0), 0,1)
      res.rand <- list(dce=tmp, dce.pvalue=tmp2)
      time.rand <- as.integer(difftime(Sys.time(), time.tmp, units="secs"))

      df.res <- data.frame(
        truth=as.vector(ground.truth$dce),
        cor=as.vector(res.cor$dce),
        pcor=as.vector(res.pcor$dce),
        dce=as.vector(res.dce$dce),
        rand=as.vector(res.rand$dce)
      )

      df.pvals <- data.frame(
        truth=as.vector(ground.truth$dce),
        cor=as.vector(res.cor$dce.pvalue),
        pcor=as.vector(res.pcor$dce.pvalue),
        dce=as.vector(res.dce$dce.pvalue),
        rand=as.vector(res.rand$dce.pvalue)
      )

      # make performance evaluation fairer by only comparing results for edges
      df.res %<>% dplyr::filter(as.vector(as(wt.graph, "matrix")) != 0)
      df.pvals %<>% dplyr::filter(as.vector(as(wt.graph, "matrix")) != 0)


      # return performance computation
      data.frame() %>%
        bind_rows(
          as.data.frame(
            cor(df.res, method="spearman", use="pairwise.complete.obs")
          ) %>%
            rownames_to_column() %>%
            dplyr::filter(rowname == "truth") %>%
            dplyr::select(-rowname, -truth) %>%
            mutate(type="correlation"),

          bind_rows(list(
            cor=get_prediction_counts(df.res$truth, df.res$cor),
            pcor=get_prediction_counts(df.res$truth, df.res$pcor),
            dce=get_prediction_counts(df.res$truth, df.res$dce),
            rand=get_prediction_counts(df.res$truth, df.res$rand)
          ), .id="name") %>%
            column_to_rownames(var="name") %>%
            t %>%
            as.data.frame %>%
            rownames_to_column(var="type"),

          data.frame(
            cor=compute.mse(df.res$truth, df.res$cor),
            pcor=compute.mse(df.res$truth, df.res$pcor),
            dce=compute.mse(df.res$truth, df.res$dce),
            rand=compute.mse(df.res$truth, df.res$rand)
          ) %>%
            mutate(type="mse"),

          data.frame(
            cor=compute.prec(df.pvals$truth, df.pvals$cor),
            pcor=compute.prec(df.pvals$truth, df.pvals$pcor),
            dce=compute.prec(df.pvals$truth, df.pvals$dce),
            rand=compute.prec(df.pvals$truth, df.pvals$rand)
          ) %>%
            mutate(type="precision"),

          data.frame(
            cor=compute.prec(df.pvals$truth, df.pvals$cor, do = "rec"),
            pcor=compute.prec(df.pvals$truth, df.pvals$pcor, do = "rec"),
            dce=compute.prec(df.pvals$truth, df.pvals$dce, do = "rec"),
            rand=compute.prec(df.pvals$truth, df.pvals$rand, do = "rec")
          ) %>%
            mutate(type="recall"),

          data.frame(
            cor=compute.prauc(df.pvals$truth, df.pvals$cor),
            pcor=compute.prauc(df.pvals$truth, df.pvals$pcor),
            dce=compute.prauc(df.pvals$truth, df.pvals$dce),
            rand=compute.prauc(df.pvals$truth, df.pvals$rand)
          ) %>%
            mutate(type="pr-auc"),

          data.frame(
            cor=compute.rocauc(df.pvals$truth, df.pvals$cor),
            pcor=compute.rocauc(df.pvals$truth, df.pvals$pcor),
            dce=compute.rocauc(df.pvals$truth, df.pvals$dce),
            rand=compute.rocauc(df.pvals$truth, df.pvals$rand)
          ) %>%
            mutate(type="roc-auc"),

          data.frame(
            cor=time.cor,
            pcor=time.pcor,
            dce=time.dce,
            rand=time.rand
          ) %>%
            mutate(type="runtime"),

          data.frame(
            cor=graph.density,
            pcor=graph.density,
            dce=graph.density,
            rand=graph.density
          ) %>%
            mutate(type="graph.density"),

          data.frame(
            cor=dispersion.estimate,
            pcor=dispersion.estimate,
            dce=dispersion.estimate,
            rand=dispersion.estimate
          ) %>%
            mutate(type="dispersion.estimate"),

          data.frame(
            cor=mean.estimate,
            pcor=mean.estimate,
            dce=mean.estimate,
            rand=mean.estimate
          ) %>%
            mutate(type="mean.estimate"),
          
        ) %>%
        mutate(parameter=parameter) %>% mutate(varied.parameter=varied.parameter)
    },
    otherwise=NULL,
    quiet=FALSE
  )
) %>% { if (file.exists("benchmark_results.csv") && append) {
            write_csv(., "benchmark_results.csv", append = append)
        } else {
            write_csv(., "benchmark_results.csv", append = append, col_names = TRUE)
        }}

if (append) {
    df.bench <- read_csv("benchmark_results.csv")
}

df.bench <- df.bench[which(df.bench$varied.parameter == varied.parameter), ]

if (!(varied.parameter %in% c("adjustment.type", "theta.fixed"))) {
    df.bench$parameter %<>% as.factor %>% fct_inseq
    df.bench %>%
    head
}

# plotting
df.bench %>%
  dplyr::filter(type == "correlation") %>%
  gather("variable", "value", -parameter, -type, -varied.parameter) %>%
ggplot(aes(x=parameter, y=value, fill=variable)) +
  geom_boxplot() +
  ylim(-1, 1) +
  ggtitle(paste("Variable:", varied.parameter)) +
  ylab("Correlation (truth vs prediction)") +
  theme_minimal(base_size=20) +
  theme(plot.title=element_text(hjust=0.5)) +
  ggsave("benchmark_correlation.pdf")

df.bench %>%
  dplyr::filter(type == "precision") %>%
  gather("variable", "value", -parameter, -type, -varied.parameter) %>%
ggplot(aes(x=parameter, y=value, fill=variable)) +
  geom_boxplot() +
  ylim(0, 1) +
  ggtitle(paste("Variable:", varied.parameter)) +
  ylab("Precision (truth vs prediction)") +
  theme_minimal(base_size=20) +
  theme(plot.title=element_text(hjust=0.5)) +
  ggsave("benchmark_precision.pdf")

df.bench %>%
  dplyr::filter(type == "recall") %>%
  gather("variable", "value", -parameter, -type, -varied.parameter) %>%
ggplot(aes(x=parameter, y=value, fill=variable)) +
  geom_boxplot() +
  ylim(0, 1) +
  ggtitle(paste("Variable:", varied.parameter)) +
  ylab("Recall (truth vs prediction)") +
  theme_minimal(base_size=20) +
  theme(plot.title=element_text(hjust=0.5)) +
  ggsave("benchmark_recall.pdf")

df.bench %>%
  dplyr::filter(type == "pr-auc") %>%
  gather("variable", "value", -parameter, -type, -varied.parameter) %>%
ggplot(aes(x=parameter, y=value, fill=variable)) +
  geom_boxplot() +
  ylim(0, 1) +
  ggtitle(paste("Variable:", varied.parameter)) +
  ylab("Precision-Recall (AUC) (truth vs prediction)") +
  theme_minimal(base_size=20) +
  theme(plot.title=element_text(hjust=0.5)) +
  ggsave("benchmark_prauc.pdf")

df.bench %>%
  dplyr::filter(type == "roc-auc") %>%
  gather("variable", "value", -parameter, -type, -varied.parameter) %>%
ggplot(aes(x=parameter, y=value, fill=variable)) +
  geom_boxplot() +
  ylim(0, 1) +
  ggtitle(paste("Variable:", varied.parameter)) +
  ylab("ROC (AUC) (truth vs prediction)") +
  theme_minimal(base_size=20) +
  theme(plot.title=element_text(hjust=0.5)) +
  ggsave("benchmark_rocauc.pdf")

df.bench %>%
  dplyr::filter(type == "mse") %>%
  gather("variable", "value", -parameter, -type, -varied.parameter) %>%
  ggplot(aes(x=parameter, y=value, fill=variable)) +
  geom_boxplot() +
  ggtitle(paste("Variable:", varied.parameter)) +
  ylab("Mean squared error") +
  theme_minimal(base_size=20) +
  theme(plot.title=element_text(hjust=0.5)) +
  ggsave("benchmark_mse.pdf")

df.bench %>%
  dplyr::filter(type == "runtime") %>%
  gather("variable", "value", -parameter, -type, -varied.parameter) %>%
  mutate(value=lubridate::as.duration(value)) %>%
ggplot(aes(x=parameter, y=value, fill=parameter)) +
  geom_boxplot() +
  scale_y_time() +
  facet_wrap(~ variable, scales="free") +
  ggtitle(paste("Variable:", varied.parameter)) +
  ylab("runtime") +
  theme_minimal() +
  theme(plot.title=element_text(hjust=0.5)) +
  ggsave("benchmark_runtime.pdf")

df.bench %>%
  dplyr::filter(grepl("^graph.", type)) %>%
  dplyr::select(dce, type, parameter, varied.parameter) %>%
ggplot(aes(x=parameter, y=dce, fill=type)) +
  geom_boxplot() +
  ggtitle(paste("Variable:", varied.parameter)) +
  ylab("value") +
  theme_minimal(base_size=20) +
  theme(plot.title=element_text(hjust=0.5)) +
  ggsave("benchmark_graph_features.pdf")

df.bench %>%
  dplyr::filter(grepl("^dispersion.", type)) %>%
  dplyr::select(dce, type, parameter, varied.parameter) %>%
ggplot(aes(x=parameter, y=dce, fill=type)) +
  geom_boxplot() +
  ggtitle(paste("Variable:", varied.parameter)) +
  ylab("value") +
  theme_minimal(base_size=20) +
  theme(plot.title=element_text(hjust=0.5)) +
  ggsave("benchmark_dispersion_estimate.pdf")

df.bench %>%
  dplyr::filter(grepl("^mean.", type)) %>%
  dplyr::select(dce, type, parameter, varied.parameter) %>%
ggplot(aes(x=parameter, y=dce, fill=type)) +
  geom_boxplot() +
  ggtitle(paste("Variable:", varied.parameter)) +
  ylab("value") +
  theme_minimal(base_size=20) +
  theme(plot.title=element_text(hjust=0.5)) +
  ggsave("benchmark_mean_estimate.pdf")
