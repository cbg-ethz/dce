library(tidyverse)

# parse commandline arguments
"
Benchmark DCE performance and runtime.

Usage:
  plotting.R
  plotting.R --input STR --output STR
  plotting.R --input STR --output STR --methods STR
  plotting.R --input STR --output STR --methods STR --parameters STR

Options:
  -h --help     Show this screen.
  --input STR   CSV file to read data from [default: benchmark_results.csv].
  --output STR  Directory to store plots in [default: plots/].
  --methods STR    Which methods to plot separated by commas, e.g., cor,pcor [default: NULL].
  --parameters STR    Which parameters to plot separated by commas, e.g., 1,5 [default: NULL].
" -> doc

arguments <- docopt::docopt(doc)


# setup environment
input.fname <- "benchmark_results.csv"
target.dir <- "plots/"

input.fname <- arguments$input
target.dir <- arguments$output
methods <- unlist(strsplit(arguments$methods,','))
parameters <- unlist(strsplit(arguments$parameters,','))

print(glue::glue("{input.fname} -> {target.dir} (methods) (parameters)"))


# helper functions
create.plots <- function(df.bench, plot.dir, varied.parameter) {
  # create performance plots
  performance.measures <- c("correlation", "mse", "precision", "recall", "f1-score", "pr-auc", "ROC-AUC","ROC-AUC (ES)")

  height <- 20
  Npar <- length(table(df.bench$parameter))
  width <- 8*Npar

  for (measure in performance.measures) {
    print(glue::glue("Plotting {measure}"))
    
    xlabel <- switch(unique(df.bench$varied.parameter),
                     "Adjustment set" = 'Confounders',
                     "Effect magnitude" = 'Maximum absolute effect size',
                     "Dispersion" = 'Dispersion strength',
                     "Latent variables" = 'Number of latent variables',
                     "Library size range" = 'Maximum library size factor',
                     "Number of samples" = 'Number of samples',
                     "Network size" = 'Number of genes in the network',
                     "Network perturbation" = 'Fraction of added/deleted edges',
                     "Prevalence of positive edges" = 'Prevelance of true differential effects',
                     "Beta distribution" = "Beta distribution",
                     "latent distribution" = "latent distribution")

    meth_order <- c('dce', 'cor', 'pcor', 'fggm', 'ldgm', 'car', 'dge', 'rand', 'dce (no library size correction)', 'dce (no latent correction)')
    meth_color <- c('red', 'lightblue', 'blue', 'orange', 'pink', 'green', 'purple', 'grey', '#ff7777', 'darkred')
    meth_color <- meth_color[meth_order %in% colnames(df.bench)]
    meth_order <- meth_order[meth_order %in% colnames(df.bench)]

    options(ggplot.discrete.fill = meth_color)

    p <- df.bench %>%
      dplyr::filter(type == measure) %>%
      gather("variable", "value", -parameter, -type, -varied.parameter, -rng.seed) %>%
      mutate(variable=factor(variable, levels=meth_order)) %>%
      ggplot(aes(x=parameter, y=value, fill=variable)) +
      geom_boxplot() +
      scale_fill_manual(values=meth_color) +
      ggtitle(paste(varied.parameter)) +
      ylab(glue::glue("{measure}")) +
      xlab(xlabel) +
      theme_minimal(base_size=20) +
      theme(plot.title=element_text(hjust=0.5)) +
      guides(fill=guide_legend(title="Methods",
                               label.them = element_text(face = 'italic')))

    if (measure %in% c("correlation")) {
      p <- p + ylim(-1, 1)
    } else if (measure %in% c("precision", "recall", "f1-score", "pr-auc", "roc-auc", "roc-auc_es", "ROC-AUC", "ROC-AUC (ES)")) {
      p <- p + ylim(0, 1)
    }

    ggsave(file.path(plot.dir, glue::glue("benchmark_{measure}.pdf")), plot = p,
           width = width, height = height, units = 'cm')
  }


  # create special plots
  p <- df.bench %>%
    dplyr::filter(type == "runtime") %>%
    gather("variable", "value", -parameter, -type, -varied.parameter, -rng.seed) %>%
    mutate(value=lubridate::as.duration(value)) %>%
    ggplot(aes(x=parameter, y=value, fill=parameter)) +
    geom_boxplot() +
    scale_y_time() +
    facet_wrap(~ variable, scales="free") +
    ggtitle(paste(varied.parameter)) +
    ylab("runtime") +
    theme_minimal() +
    theme(plot.title=element_text(hjust=0.5))
    ggsave(file.path(plot.dir, "benchmark_runtime.pdf"),
           width = width, height = height, units = 'cm', plot = p)

  meth <- colnames(df.bench)[1]

  p <- df.bench %>%
    dplyr::filter(grepl("^dce.", type)) %>%
    dplyr::select(!!colnames(df.bench)[1], type, parameter, varied.parameter, meth) %>%
    ggplot(aes_string(x='parameter', y=meth, fill='type')) +
    scale_y_continuous(trans = 'log10') +
    # scale_y_log10() +
    geom_boxplot() +
    ggtitle(paste(varied.parameter)) +
    ylab("value") +
    theme_minimal(base_size=20) +
    theme(plot.title=element_text(hjust=0.5))
    ggsave(file.path(plot.dir, "benchmark_dce_range.pdf"),
           width = width, height = height, units = 'cm', plot = p)

  p <- df.bench %>%
    dplyr::filter(grepl("^graph.", type)) %>%
    dplyr::select(!!colnames(df.bench)[1], type, parameter, varied.parameter) %>%
    ggplot(aes_string(x='parameter', y=meth, fill='type')) +
    geom_boxplot() +
    ggtitle(paste(varied.parameter)) +
    ylab("value") +
    theme_minimal(base_size=20) +
    theme(plot.title=element_text(hjust=0.5))
    ggsave(file.path(plot.dir, "benchmark_graph_features.pdf"),
           width = width, height = height, units = 'cm', plot = p)

  p <- df.bench %>%
    dplyr::filter(grepl("^lib.", type)) %>%
    dplyr::select(!!colnames(df.bench)[1], type, parameter, varied.parameter) %>%
    ggplot(aes_string(x='parameter', y=meth, fill='type')) +
    geom_boxplot() +
    ggtitle(paste(varied.parameter)) +
    ylab("value") +
    theme_minimal(base_size=20) +
    theme(plot.title=element_text(hjust=0.5))
    ggsave(file.path(plot.dir, "benchmark_lib_size_stats.pdf"),
           width = width, height = height, units = 'cm', plot = p)

  p <- df.bench %>%
    dplyr::filter(grepl("^dispersion.", type)) %>%
    dplyr::select(!!colnames(df.bench)[1], type, parameter, varied.parameter) %>%
    ggplot(aes_string(x='parameter', y=meth, fill='type')) +
    geom_boxplot() +
    ggtitle(paste(varied.parameter)) +
    ylab("value") +
    theme_minimal(base_size=20) +
    theme(plot.title=element_text(hjust=0.5))
    ggsave(file.path(plot.dir, "benchmark_dispersion_estimate.pdf"),
           width = width, height = height, units = 'cm', plot = p)

  p <- df.bench %>%
    dplyr::filter(grepl("^mean.", type)) %>%
    dplyr::select(!!colnames(df.bench)[1], type, parameter, varied.parameter) %>%
    ggplot(aes_string(x='parameter', y=meth, fill='type')) +
    geom_boxplot() +
    ggtitle(paste(varied.parameter)) +
    ylab("value") +
    theme_minimal(base_size=20) +
    theme(plot.title=element_text(hjust=0.5))
    ggsave(file.path(plot.dir, "benchmark_mean_estimate.pdf"),
           width = width, height = height, units = 'cm', plot = p)

  p <- df.bench %>%
    dplyr::filter(grepl("^prevalence$", type)) %>%
    dplyr::select(!!colnames(df.bench)[1], type, parameter, varied.parameter) %>%
    ggplot(aes_string(x='parameter', y=meth, fill='type')) +
    geom_boxplot() +
    ggtitle(paste(varied.parameter)) +
    ylab("value") +
    theme_minimal(base_size=20) +
    theme(plot.title=element_text(hjust=0.5))
    ggsave(file.path(plot.dir, "benchmark_prevalence.pdf"),
           width = width, height = height, units = 'cm', plot = p)
}


# read data
if (!dir.exists(target.dir)) {
  dir.create(target.dir, recursive = TRUE)
}

df.bench <- read_csv(input.fname)

ld.idx <- which(df.bench$varied.parameter == "latent.dist")
df.bench$parameter[ld.idx] <- mgsub::mgsub(df.bench$parameter[ld.idx], c("1","2"), c("uniform", "exponential"))

df.bench <- df.bench[!is.na(df.bench$varied.parameter), ]
if (methods[1]!='NULL') {
  df.bench <- df.bench[,colnames(df.bench) %in% c(methods,'type','parameter','varied.parameter','rng.seed')]
}
df.bench <- df.bench[, apply(df.bench[df.bench$type == 'roc-auc', ], 2, function(x) {
  return(!all(is.na(x)))
})]

if (parameters[1]!='NULL') {
  df.bench <- df.bench[df.bench$parameter %in% parameters,]
}

# rename variables for paper ready (;)) figures:
df.bench$varied.parameter <- mgsub::mgsub(df.bench$varied.parameter,
                                          c("adjustment.type", "beta.magnitude", "dispersion", "latent", "lib.size.range", "mt.samples", "node.num", "perturb", "true.positives","beta.dist","latent.dist"),
                                          c("Adjustment set", "Effect magnitude", "Dispersion", "Latent variables", "Library size range", "Number of samples", "Network size", "Network perturbation", "Prevalence of positive edges", "Beta distribution",
                                            "latent distribution"))
colnames(df.bench) <- mgsub::mgsub(colnames(df.bench),
                                   c('dce.lm.tpm', 'fggm', 'cor', 'pcorz', 'dce.nolib', 'dce.lm.tpm.nolatent'),
                                   c('dce', 'fggm', 'cor', 'pcor', 'dce (no library size correction)', 'dce (no latent correction)'))
if (!'dce' %in% colnames(df.bench)) {
  colnames(df.bench) <- gsub('.HC', '', colnames(df.bench))
}
df.bench$type <- mgsub::mgsub(df.bench$type,
                              c('roc-auc','roc-auc_es'),
                              c('ROC-AUC','ROC-AUC (ES)'))

tmp <- df.bench %>% pull(varied.parameter) %>% unique
if (length(tmp) != 1) {
  tmp <- table(df.bench$varied.parameter)
  varied.parameter <- names(tmp)[which.max(tmp)]
  df.bench <- df.bench[!is.na(df.bench$parameter),]
} else {
  varied.parameter <- tmp
}
if (varied.parameter == 'Dispersion') {
  df.bench$parameter <- 1/df.bench$parameter
}

if (!any(is.na(as.numeric(df.bench$parameter)))) {
  df.bench$parameter %<>% as.factor %>% fct_inseq
}

# create plots
create.plots(df.bench, target.dir, varied.parameter)
