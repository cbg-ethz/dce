create.plots <- function(df.bench, target.dir) {
  # create performance plots
  performance.measures <- c("correlation", "mse", "precision", "recall", "pr-auc", "roc-auc")

  for (measure in performance.measures) {
    print(glue::glue("Plotting {measure}"))

    df.bench %>%
      dplyr::filter(type == measure) %>%
      gather("variable", "value", -parameter, -type, -varied.parameter) %>%
      ggplot(aes(x=parameter, y=value, fill=variable)) +
      geom_boxplot() +
      ylim(-1, 1) +
      ggtitle(paste("Variable:", varied.parameter)) +
      ylab(glue::glue("{measure} (truth vs prediction)")) +
      theme_minimal(base_size=20) +
      theme(plot.title=element_text(hjust=0.5)) +
      ggsave(file.path(target.dir, glue::glue("benchmark_{measure}.pdf")))
  }


  # create special plots
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
    ggsave(file.path(target.dir, "benchmark_runtime.pdf"))

  df.bench %>%
    dplyr::filter(grepl("^graph.", type)) %>%
    dplyr::select(dce, type, parameter, varied.parameter) %>%
    ggplot(aes(x=parameter, y=dce, fill=type)) +
    geom_boxplot() +
    ggtitle(paste("Variable:", varied.parameter)) +
    ylab("value") +
    theme_minimal(base_size=20) +
    theme(plot.title=element_text(hjust=0.5)) +
    ggsave(file.path(target.dir, "benchmark_graph_features.pdf"))

  df.bench %>%
    dplyr::filter(grepl("^dispersion.", type)) %>%
    dplyr::select(dce, type, parameter, varied.parameter) %>%
    ggplot(aes(x=parameter, y=dce, fill=type)) +
    geom_boxplot() +
    ggtitle(paste("Variable:", varied.parameter)) +
    ylab("value") +
    theme_minimal(base_size=20) +
    theme(plot.title=element_text(hjust=0.5)) +
    ggsave(file.path(target.dir, "benchmark_dispersion_estimate.pdf"))

  df.bench %>%
    dplyr::filter(grepl("^mean.", type)) %>%
    dplyr::select(dce, type, parameter, varied.parameter) %>%
    ggplot(aes(x=parameter, y=dce, fill=type)) +
    geom_boxplot() +
    ggtitle(paste("Variable:", varied.parameter)) +
    ylab("value") +
    theme_minimal(base_size=20) +
    theme(plot.title=element_text(hjust=0.5)) +
    ggsave(file.path(target.dir, "benchmark_mean_estimate.pdf"))
}
