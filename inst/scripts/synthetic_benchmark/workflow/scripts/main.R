library(tidyverse)
library(magrittr)
library(graph)
library(naturalsort)

devtools::load_all("../../")

source("workflow/scripts/helper_functions.R")
source("workflow/scripts/models.R")
source("workflow/scripts/performance_measures.R")
source("workflow/scripts/LDGM.R")
source("workflow/scripts/FastGGM.R")
source("workflow/scripts/Carnival.R")
source("workflow/scripts/DGE.R")

# parse commandline arguments
"
Benchmark DCE performance and runtime.

Usage:
  main.R
  main.R --variable NAME --values VALUES
  main.R --variable NAME --values VALUES --methods STR
  main.R --variable NAME --values VALUES --replicates INT --methods STR
  main.R --variable NAME --values VALUES --replicates INT --output STR
  main.R --variable NAME --values VALUES --replicates INT --output STR --methods STR
  main.R --variable NAME --values VALUES --append BOOL --replicates INT --link STR --output STR

Options:
  -h --help        Show this screen.
  --variable NAME  Which property to vary [default: node.num].
  --values VALUES  What values to assign to varied property [default: 20,50,100].
  --append BOOL    If TRUE appends the results of this/these run(s) to an existing results file [default: FALSE].
  --replicates INT Number of simulation runs [default: 100].
  --link STR       Either log or identity as link function [default: identity].
  --output STR     CSV File to store results in [default: benchmark_results.csv].
  --methods STR    Which methods to benchmark [default: NULL].
" -> doc

arguments <- docopt::docopt(doc)


# global parameters
node.num <- 100
wt.samples <- 200
mt.samples <- 200

beta.magnitude <- 1
beta.dist <- 1
latent.dist <- 1
dist.mean <- 100
dispersion <- 1
adjustment.type <- "parents"
effect.type <- "total"

sample.kegg <- FALSE
append <- FALSE

perturb <- 0
true.positives <- 0.5
lib.size.range <- 10
latent <- 0


# special parameters which can later be modified from commandline
output.fname <- "benchmark_results.csv"
replicate.count <- 100
link.method <- "identity"
methods <- NULL

# parse parameters
varied.parameter <- arguments$variable
parameter.list <- unlist(
  purrr::map(strsplit(arguments$values, ",")[[1]], type.convert)
)

if (varied.parameter == "latent.dist") {
  latent <- 10
}

output.fname <- arguments$output
replicate.count <- as.numeric(arguments$replicates)
append <- as.logical(arguments$append)

if (arguments$methods != "NULL") {
  methods <- strsplit(arguments$methods, ",")[[1]]
}

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
  kegg.dag <- readRDS("resources/pathways.rds")
  node.num <- 10^9
  replicate.count <- length(kegg.dag)
}


# run benchmark
df.bench <- purrr::pmap_dfr(
  list(parameter=rep(parameter.list, each=replicate.count), index=rep(seq_len(replicate.count), length(parameter.list))),
  purrr::possibly(
    function(parameter, index) {
      # handle parameterization
      rng.seed <- seed.list[index]
      set.seed(rng.seed)

      switch(
        varied.parameter,
        node.num = { node.num <- parameter },
        wt.samples = { wt.samples <- parameter },
        mt.samples = { mt.samples <- parameter },
        beta.magnitude = { beta.magnitude <- parameter },
        beta.dist = { beta.dist <- parameter },
        latent.dist = { latent.dist <- parameter },
        dispersion = { dispersion <- parameter },
        adjustment.type = { adjustment.type <- parameter },
        effect.type = { effect.type <- parameter },
        perturb = { perturb <- parameter },
        true.positives = { true.positives <- parameter },
        lib.size.range = { lib.size.range <- parameter },
        latent = { latent <- parameter },

        total.samples={
          wt.samples <- round(parameter / 2)
          mt.samples <- round(parameter / 2)
        }
      )

      print(glue::glue("seed={rng.seed} node.num={node.num} wt.samples={wt.samples} mt.samples={mt.samples} beta.magnitude={beta.magnitude} beta.dist={beta.dist} latent.dist={latent.dist} dispersion={dispersion} adjustment.type={adjustment.type} effect.type={effect.type} perturb={perturb} true.positives={true.positives} lib.size.range={lib.size.range} latent={latent}"))


      # generate graphs
      if (sample.kegg) {
        graphs <- sample.graph.from.kegg(kegg.dag)
      } else {
        if (beta.dist == 1) {
          betameth <- "unif"
        } else {
          betameth <- "exp"
        }
        graphs <- generate.random.graphs(node.num, beta.magnitude,
                                         true.positives, max_par = 10,
                                         mineff = 0, method = betameth)
      }

      wt.graph <- graphs$wt
      mt.graph <- graphs$mt

      prevalence <- compute.prevalence(wt.graph, mt.graph)

      # compute dce stats
      dce.stats <- compute.dce.stats(wt.graph, mt.graph)

      # generate data
      pop.size <- 10000
      if (latent.dist == 1) {
        latent.fun <- "unif"
      } else {
        latent.fun <- "exp"
      }
      wt.X <- simulate_data(wt.graph, n = wt.samples, dist_dispersion = dispersion, dist_mean = dist.mean, pop_size = pop.size, latent = latent, latent.fun = latent.fun)
      mt.X <- simulate_data(mt.graph, n = mt.samples, dist_dispersion = dispersion, dist_mean = dist.mean, pop_size = pop.size, latent = latent, latent.fun = latent.fun)

      # library size difference
      lib.size.mean <- (lib.size.range+1)/2
      lib.size.sd <- lib.size.range/10
      ptruncnorm <- dnorm(1:lib.size.range,lib.size.mean,lib.size.sd)/(pnorm(lib.size.range+1,lib.size.mean,lib.size.sd)-pnorm(0,lib.size.mean,lib.size.sd))
      lib.size.gtn <- sample(1:lib.size.range,nrow(wt.X)+nrow(mt.X),replace=TRUE,prob=ptruncnorm)
      wt.X <- wt.X*lib.size.gtn[seq_len(wt.samples)]
      mt.X <- mt.X*lib.size.gtn[(wt.samples+1):(wt.samples+mt.samples)]

      xt <- c(rep(0,nrow(wt.X)),rep(1,nrow(mt.X)))
      names(xt) <- "group"
      design <- model.matrix(~xt)
      dispersion.estimate <- NA # estimateTheta(rbind(wt.X, mt.X), design = design)
      mean.estimate <- mean(rbind(wt.X, mt.X))

      # check how close we get to library size
      lib.size <- apply(rbind(wt.X, mt.X), 1, sum)
      lib.size <- round(lib.size/(10^min(round(log10(lib.size)))))
      lib.size.stats <- cor(lib.size, lib.size.gtn)


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
      latent2 <- FALSE
      if (varied.parameter == 'latent' | varied.parameter == 'latent.dist') {
        latent2 <- 'auto'
      }
      res <- run.all.models(
        wt.graph, wt.X,
        mt.graph, mt.X,
        wt.graph.perturbed,
        beta.magnitude,
        methods = methods,
        effect.type = effect.type,
        adjustment.type = adjustment.type,
        latent = latent2
      )

      if (is.null(methods)) {
        methods <- colnames(res$runtime)
      }

      df.edges <- res$edges
      df.pvalues <- res$pvalues
      df.runtime <- res$runtime

      # modify predictions
      #  * if edge exists only in original graph (but not in perturbed one), it should be a false negative if truth != 0
      #  * for performance evaluation use all entries with edge in original or perturbed graph
      tmp.graph <- as_adjmat(wt.graph)
      tmp.graph.perturbed <- as_adjmat(wt.graph.perturbed)
      df.all <- bind_cols(
        data.frame(orig.edge = as.numeric(as.vector(tmp.graph) != 0)),
        data.frame(pert.edge = as.numeric(as.vector(tmp.graph.perturbed) != 0)),
        df.pvalues
      )
      df.all.edges <- bind_cols(
        data.frame(orig.edge = as.numeric(as.vector(tmp.graph) != 0)),
        data.frame(pert.edge = as.numeric(as.vector(tmp.graph.perturbed) != 0)),
        df.edges
      )

      for (method in methods) {
        df.all[[method]] <- case_when(
          df.all$pert.edge == df.all$orig.edge | df.all$pert.edge == 1 ~ df.all[[method]],
          df.all$pert.edge != df.all$orig.edge ~ 1
        )
      }

      df.all %<>%
        dplyr::filter(orig.edge | pert.edge)
      df.all.edges %<>%
        dplyr::filter(orig.edge | pert.edge)

      df.pvalues.mod <- df.all %>%
        select(-orig.edge, -pert.edge)
      df.edges.mod <- df.all.edges %>%
        select(-orig.edge, -pert.edge)

      apply.performance.measure(df.edges.mod, methods, compute.rocauc_es, "roc-auc_es")
      apply.performance.measure(df.pvalues.mod, methods, compute.prauc, "pr-auc")
      apply.performance.measure(df.pvalues.mod, methods, compute.rocauc, "roc-auc")
      cor(df.edges[c("truth", methods)], method = "spearman", use = "pairwise.complete.obs")

      # return performance computation
      data.frame() %>%
        bind_rows(
          as.data.frame(
            cor(df.edges.mod[c("truth", methods)], method = "spearman", use = "pairwise.complete.obs")
          ) %>%
            rownames_to_column() %>%
            dplyr::filter(rowname == "truth") %>%
            dplyr::select(-rowname, -truth) %>%
            mutate(type = "correlation"),

          purrr::map_dfr(methods, function(method) {
            get.classification.counts(df.pvalues.mod, method) %>%
              as.data.frame %>%
              mutate(name = method)
          }) %>%
            column_to_rownames(var = "name") %>%
            t %>%
            as.data.frame %>%
            rownames_to_column(var = "type"),

          apply.performance.measure(df.edges.mod, methods, compute.mse, "mse"),
          apply.performance.measure(df.pvalues.mod, methods, compute.precision, "precision"),
          apply.performance.measure(df.pvalues.mod, methods, compute.recall, "recall"),
          apply.performance.measure(df.pvalues.mod, methods, compute.f1score, "f1-score"),
          apply.performance.measure(df.pvalues.mod, methods, compute.prauc, "pr-auc"),
          apply.performance.measure(df.pvalues.mod, methods, compute.rocauc, "roc-auc"),
          apply.performance.measure(df.edges.mod, methods, compute.rocauc_es, "roc-auc_es"),

          df.runtime %>% mutate(type = "runtime"),

          methods %>% purrr::map_dfc(setNames, object = list(graph.density)) %>% mutate(type = "graph.density"),
          methods %>% purrr::map_dfc(setNames, object = list(lib.size.stats)) %>% mutate(type = "lib.size.stats"),
          methods %>% purrr::map_dfc(setNames, object = list(dce.stats$min)) %>% mutate(type = "dce.min"),
          methods %>% purrr::map_dfc(setNames, object = list(dce.stats$max)) %>% mutate(type = "dce.max"),
          methods %>% purrr::map_dfc(setNames, object = list(dce.stats$median)) %>% mutate(type = "dce.median"),
          methods %>% purrr::map_dfc(setNames, object = list(dce.stats$mean)) %>% mutate(type = "dce.mean"),
          methods %>% purrr::map_dfc(setNames, object = list(dispersion.estimate)) %>% mutate(type = "dispersion.estimate"),
          methods %>% purrr::map_dfc(setNames, object = list(mean.estimate)) %>% mutate(type = "mean.estimate"),
          methods %>% purrr::map_dfc(setNames, object = list(prevalence)) %>% mutate(type = "prevalence")
        ) %>%
          mutate(parameter = parameter, varied.parameter = varied.parameter, rng.seed = rng.seed)
    },
    otherwise = NULL,
    quiet = FALSE
  )
) %>%
  write_csv(output.fname, append = file.exists(output.fname))
