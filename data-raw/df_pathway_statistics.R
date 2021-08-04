devtools::load_all()

df_pathway_statistics <- get_pathway_info(include_network_statistics = TRUE)

usethis::use_data(df_pathway_statistics, overwrite = TRUE)
