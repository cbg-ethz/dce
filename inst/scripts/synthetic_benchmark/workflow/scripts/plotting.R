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
  performance.measures <- c("correlation", "mse", "precision", "recall", "f1-score", "pr-auc", "roc-auc")

  for (measure in performance.measures) {
    print(glue::glue("Plotting {measure}"))

    p <- df.bench %>%
      dplyr::filter(type == measure) %>%
      gather("variable", "value", -parameter, -type, -varied.parameter, -rng.seed) %>%
      ggplot(aes(x=parameter, y=value, fill=variable)) +
      geom_boxplot() +
      ggtitle(paste("Variable:", varied.parameter)) +
      ylab(glue::glue("{measure}")) +
      theme_minimal(base_size=20) +
      theme(plot.title=element_text(hjust=0.5))

    if (measure %in% c("correlation")) {
      p <- p + ylim(-1, 1)
    } else if (measure %in% c("precision", "recall", "f1-score", "pr-auc", "roc-auc")) {
      p <- p + ylim(0, 1)
    }

    ggsave(file.path(plot.dir, glue::glue("benchmark_{measure}.pdf")), plot = p)
  }


  # create special plots
  df.bench %>%
    dplyr::filter(type == "runtime") %>%
    gather("variable", "value", -parameter, -type, -varied.parameter, -rng.seed) %>%
    mutate(value=lubridate::as.duration(value)) %>%
    ggplot(aes(x=parameter, y=value, fill=parameter)) +
    geom_boxplot() +
    scale_y_time() +
    facet_wrap(~ variable, scales="free") +
    ggtitle(paste("Variable:", varied.parameter)) +
    ylab("runtime") +
    theme_minimal() +
    theme(plot.title=element_text(hjust=0.5)) +
    ggsave(file.path(plot.dir, "benchmark_runtime.pdf"))

  df.bench %>%
    dplyr::filter(grepl("^dce.", type)) %>%
    dplyr::select(cor, type, parameter, varied.parameter, cor) %>%
    ggplot(aes(x=parameter, y=cor, fill=type)) +
    scale_y_continuous(trans = 'log10') +
    # scale_y_log10() +
    geom_boxplot() +
    ggtitle(paste("Variable:", varied.parameter)) +
    ylab("value") +
    theme_minimal(base_size=20) +
    theme(plot.title=element_text(hjust=0.5)) +
    ggsave(file.path(plot.dir, "benchmark_dce_range.pdf"))

  df.bench %>%
    dplyr::filter(grepl("^graph.", type)) %>%
    dplyr::select(cor, type, parameter, varied.parameter) %>%
    ggplot(aes(x=parameter, y=cor, fill=type)) +
    geom_boxplot() +
    ggtitle(paste("Variable:", varied.parameter)) +
    ylab("value") +
    theme_minimal(base_size=20) +
    theme(plot.title=element_text(hjust=0.5)) +
    ggsave(file.path(plot.dir, "benchmark_graph_features.pdf"))

  df.bench %>%
    dplyr::filter(grepl("^lib.", type)) %>%
    dplyr::select(cor, type, parameter, varied.parameter) %>%
    ggplot(aes(x=parameter, y=cor, fill=type)) +
    geom_boxplot() +
    ggtitle(paste("Variable:", varied.parameter)) +
    ylab("value") +
    theme_minimal(base_size=20) +
    theme(plot.title=element_text(hjust=0.5)) +
    ggsave(file.path(plot.dir, "benchmark_lib_size_stats.pdf"))

  df.bench %>%
    dplyr::filter(grepl("^dispersion.", type)) %>%
    dplyr::select(cor, type, parameter, varied.parameter) %>%
    ggplot(aes(x=parameter, y=cor, fill=type)) +
    geom_boxplot() +
    ggtitle(paste("Variable:", varied.parameter)) +
    ylab("value") +
    theme_minimal(base_size=20) +
    theme(plot.title=element_text(hjust=0.5)) +
    ggsave(file.path(plot.dir, "benchmark_dispersion_estimate.pdf"))

  df.bench %>%
    dplyr::filter(grepl("^mean.", type)) %>%
    dplyr::select(cor, type, parameter, varied.parameter) %>%
    ggplot(aes(x=parameter, y=cor, fill=type)) +
    geom_boxplot() +
    ggtitle(paste("Variable:", varied.parameter)) +
    ylab("value") +
    theme_minimal(base_size=20) +
    theme(plot.title=element_text(hjust=0.5)) +
    ggsave(file.path(plot.dir, "benchmark_mean_estimate.pdf"))

  df.bench %>%
    dplyr::filter(grepl("^prevalence$", type)) %>%
    dplyr::select(cor, type, parameter, varied.parameter) %>%
    ggplot(aes(x=parameter, y=cor, fill=type)) +
    geom_boxplot() +
    ggtitle(paste("Variable:", varied.parameter)) +
    ylab("value") +
    theme_minimal(base_size=20) +
    theme(plot.title=element_text(hjust=0.5)) +
    ggsave(file.path(plot.dir, "benchmark_prevalence.pdf"))
}


# read data
if (!dir.exists(target.dir)) {
  dir.create(target.dir, recursive = TRUE)
}

df.bench <- read_csv(input.fname)
if (methods[1]!='NULL') {
    df.bench <- df.bench[,colnames(df.bench) %in% c(methods,'type','parameter','varied.parameter','rng.seed')]
}

if (parameters[1]!='NULL') {
  df.bench <- df.bench[df.bench$parameter %in% parameters,]
}

# check amount of runs and NAs
print('runs:')
print(table(df.bench$parameter[df.bench$type=='correlation']))
print('dce NAs for correlation:')
print(summary(df.bench$dce.lm.tpm[df.bench$type=='correlation']))
print(summary(df.bench$dce.tpm[df.bench$type=='correlation']))
print('dce NAs for roc-auc:')
print(summary(df.bench$dce.lm.tpm[df.bench$type=='roc-auc']))
print(summary(df.bench$dce.tpm[df.bench$type=='roc-auc']))

tmp <- df.bench %>% pull(varied.parameter) %>% unique
if (length(tmp) != 1) {
  tmp <- table(df.bench$varied.parameter)
  varied.parameter <- names(tmp)[which.max(tmp)]
  df.bench <- df.bench[!is.na(df.bench$parameter),]
} else {
  varied.parameter <- tmp
}

if (!any(is.na(as.numeric(df.bench$parameter)))) {
  df.bench$parameter %<>% as.factor %>% fct_inseq
}


# create plots
create.plots(df.bench, target.dir, varied.parameter)
