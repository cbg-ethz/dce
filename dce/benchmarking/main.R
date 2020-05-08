library(tidyverse)
library(magrittr)
library(graph)

devtools::load_all("..")

source("helper_functions.R")
source("models.R")
source("performance_measures.R")


# parse commandline arguments
"
Benchmark DCE performance and runtime.

Usage:
  main.R
  main.R --variable NAME --values VALUES --append BOOL --replicates INT --link STR --output STR

Options:
  -h --help        Show this screen.
  --variable NAME  Which property to vary [default: node.num].
  --values VALUES  What values to assign to varied property [default: 20,50,100].
  --append BOOL    If TRUE appends the results of this/these run(s) to an existing results file [default: FALSE].
  --replicates INT Number of simulation runs [default: 100].
  --link STR       Either log or identity as link function [default: identity].
  --output STR     CSV File to store results in [default: benchmark_results.csv].
" -> doc

arguments <- docopt::docopt(doc)


# global parameters
node.num <- 100
wt.samples <- 200
mt.samples <- 200

beta.magnitude <- 1
dist.mean <- 1000
dispersion <- 1
adjustment.type <- "parents"

sample.kegg <- FALSE
append <- FALSE

perturb <- 0
true.positives <- 0.5


# special parameters which can later be modified from commandline
output.fname <- "benchmark_results.csv"
replicate.count <- 100
link.method <- "identity"


# parse parameters
varied.parameter <- arguments$variable
parameter.list <- unlist(
  purrr::map(strsplit(arguments$values, ",")[[1]], type.convert)
)

output.fname <- arguments$output
replicate.count <- as.numeric(arguments$replicates)
append <- as.logical(arguments$append)

link.method <- arguments$link
if (link.method == "log") {
  beta.magnitude <- beta.magnitude * 0.001
}

print(glue::glue("Benchmark parameters:"))
print(glue::glue("  Varied parameter: {varied.parameter}"))
print(glue::glue("  Parameter: {parameter.list}"))


# further preparations
seed.list <- sample(seq_len(10^9), replicate.count)

if (sample.kegg) {
  kegg.dag <- readRDS("pathways.rds")
  node.num <- 10^9
  replicate.count <- length(kegg.dag)
}

if (link.method == "log") {
  link <- function(x, offset = 0) { exp(x + offset) }
} else {
  link <- negative.binomial.special()$linkfun
}


# run benchmark
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

      set.seed(seed.list[index])
      print(glue::glue("seed={seed.list[index]} node.num={node.num} wt.samples={wt.samples} mt.samples={mt.samples} beta.magnitude={beta.magnitude} dispersion={dispersion} adjustment.type={adjustment.type} perturb={perturb}"))


      # generate graphs
      if (sample.kegg) {
        graphs <- sample.graph.from.kegg(kegg.dag)
      } else {
          graphs <- generate.random.graphs(node.num, beta.magnitude, true.positives)
      }

      wt.graph <- graphs$wt
      mt.graph <- graphs$mt


      # generate data
      wt.X <- simulate_data(wt.graph, n = wt.samples, dist.dispersion = dispersion, dist.mean = dist.mean)
      mt.X <- simulate_data(mt.graph, n = mt.samples, dist.dispersion = dispersion, dist.mean = dist.mean)

      dispersion.estimate <- estimateTheta(rbind(wt.X, mt.X))
      mean.estimate <- mean(rbind(wt.X, mt.X))


      # sanity checks
      if (any(is.nan(wt.X)) || any(is.nan(mt.X))) {
        stop("Malformed simulated data")
      }


      # perturb dag
      wt.graph.perturbed <- perturb.dag(wt.graph, perturb)


      # compute graph features
      tmp <- as(wt.graph, "matrix")
      tmp[which(tmp != 0)] <- 1
      graph.density <- sum(tmp) / ((dim(tmp)[1] * (dim(tmp)[1] - 1)) / 2)


      # run models
      res <- run.all.models(
        wt.graph, wt.X,
        mt.graph, mt.X,
        wt.graph.perturbed,
        beta.magnitude
      )

      df.edges <- res$edges
      df.pvalues <- res$pvalues
      df.runtime <- res$runtime


      # make performance evaluation fairer by only comparing results for existing edges
      df.edges %<>% dplyr::filter(as.vector(as(wt.graph, "matrix")) != 0)
      df.pvalues %<>% dplyr::filter(as.vector(as(wt.graph, "matrix")) != 0)


      # return performance computation
      data.frame() %>%
        bind_rows(
          as.data.frame(
            cor(df.edges, method="spearman", use="pairwise.complete.obs")
          ) %>%
            rownames_to_column() %>%
            dplyr::filter(rowname == "truth") %>%
            dplyr::select(-rowname, -truth) %>%
            mutate(type="correlation"),

          bind_rows(list(
            cor=get_prediction_counts(df.edges$truth, df.edges$cor),
            pcor=get_prediction_counts(df.edges$truth, df.edges$pcor),
            dce=get_prediction_counts(df.edges$truth, df.edges$dce),
            dce.lr=get_prediction_counts(df.edges$truth, df.edges$dce.lr),
            rand=get_prediction_counts(df.edges$truth, df.edges$rand)
          ), .id="name") %>%
            column_to_rownames(var="name") %>%
            t %>%
            as.data.frame %>%
            rownames_to_column(var="type"),

          apply.performance.measure(df.edges, compute.mse, "mse"),
          apply.performance.measure(df.pvalues, compute.prec, "precision"),
          apply.performance.measure(df.pvalues, compute.prec, "recall", do = "rec"),
          apply.performance.measure(df.pvalues, compute.prauc, "pr-auc"),
          apply.performance.measure(df.pvalues, compute.rocauc, "roc-auc"),

          df.runtime %>% mutate(type="runtime"),

          data.frame(
            cor=graph.density,
            pcor=graph.density,
            dce=graph.density,
            dce.lr=graph.density,
            rand=graph.density
          ) %>%
            mutate(type="graph.density"),

          data.frame(
            cor=dispersion.estimate,
            pcor=dispersion.estimate,
            dce=dispersion.estimate,
            dce.lr=dispersion.estimate,
            rand=dispersion.estimate
          ) %>%
            mutate(type="dispersion.estimate"),

          data.frame(
            cor=mean.estimate,
            pcor=mean.estimate,
            dce=mean.estimate,
            dce.lr=mean.estimate,
            rand=mean.estimate
          ) %>%
            mutate(type="mean.estimate"),

        ) %>%
          mutate(parameter=parameter) %>% mutate(varied.parameter=varied.parameter)
    },
    otherwise = NULL,
    quiet = FALSE
  )
) %>% {
  if (append && file.exists(output.fname)) {
    write_csv(., output.fname, append = append)
  } else {
    write_csv(., output.fname, append = append, col_names = TRUE)
  }
}
