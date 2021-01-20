test_that("general idea works", {
  set.seed(42)

  # create pathway
  graph <- as(dce::create_random_DAG(node_num = 20, prob = .2), "matrix")

  # generate data
  X <- simulate_data(graph, n = 10000, latent = 0)

  latent_count <- 4
  X_confounded <- simulate_data(graph, n = 10000, latent = latent_count)

  # parent adjustment set
  node_name <- "n6"
  parent_set <- names(which(graph[, node_name] != 0))

  expect_gt(length(parent_set), 0)

  # extract data
  df <- data.frame(X[, c(node_name, parent_set)])
  df_confounded <- data.frame(X_confounded[, c(node_name, parent_set)])

  # account for confounding
  fit_pca <- prcomp(scale(X_confounded))
  plot(fit_pca$sdev)
  fit_scores <- fit_pca$x[, 1:latent_count]

  df_confounded_pca <- cbind(df_confounded, fit_scores)
  df_confounded_pca %>% head

  # fit models
  fit <- lm(paste0(node_name, " ~ ."), data = df)
  fit_confounded <- lm(paste0(node_name, " ~ ."), data = df_confounded)
  fit_confounded_pca <- lm(paste0(node_name, " ~ ."), data = df_confounded_pca)

  # evaluate results
  true_value <- graph[parent_set, node_name]
  names(true_value) <- parent_set

  expect_equal(true_value, coef(fit)[parent_set], tolerance = 1e-2)
  expect_false(
    isTRUE(
      all.equal(true_value, coef(fit_confounded)[parent_set], tolerance = 1e-2)
    )
  )

  # TODO: test `coef(fit_confounded_pca)[parent_set]`
})


test_that("deconfounding can improve fit", {
  set.seed(42)

  graph_wt <- dce::create_random_DAG(node_num = 10, prob = .8)
  X_wt <- simulate_data(graph_wt, n = 10000, latent = 10)

  graph_mt <- dce::resample_edge_weights(graph_wt)
  X_mt <- simulate_data(graph_mt, n = 10000, latent = 10)

  res <- dce::dce(graph_wt, X_wt, X_mt, solver = "lm")

  res %>%
    as.data.frame %>%
    mutate(
      truth = as.vector(
        dce::trueEffects(graph_mt)
      ) - as.vector(
        dce::trueEffects(graph_wt)
      )
    ) %>%
    dplyr::select(dce, truth)

  # TODO: test result
})
