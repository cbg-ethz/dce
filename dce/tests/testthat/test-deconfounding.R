test_that("general idea works", {
  set.seed(42)

  # create pathway
  graph <- as(dce::create_random_DAG(n = 20, prob = .2), "matrix")
  # dce::plot_network(graph)

  # generate data
  X <- simulate_data(graph, n = 10000, latent = 0)

  latent_count <- 4
  X_confounded <- simulate_data(graph, n = 10000, latent = latent_count)

  # parent adjustment set
  node_name = "n6"
  parent_set = names(which(graph[, node_name] != 0))

  expect_gt(length(parent_set), 0)

  # extract data
  df <- data.frame(X[, c(node_name, parent_set)])
  df_confounded <- data.frame(X_confounded[, c(node_name, parent_set)])

  # account for confounding
  fit_pca = prcomp(scale(X_confounded))
  # plot(fit_pca)
  # plot(fit_pca$sdev)
  fit_scores = fit_pca$x[, 1:latent_count]

  df_confounded_pca = cbind(df_confounded, fit_scores)

  # fit models
  fit <- lm(paste0(node_name ,' ~ .'), data = df)
  fit_confounded <- lm(paste0(node_name ,' ~ .'), data = df_confounded)
  fit_confounded_pca <- lm(paste0(node_name ,'~ .'), data = df_confounded_pca)

  # evaluate results
  graph[parent_set, node_name]
  coef(fit)[-1]
  coef(fit_confounded)[-1]
  coef(fit_confounded_pca)[-1] # c(-1, -6:-9)
})
