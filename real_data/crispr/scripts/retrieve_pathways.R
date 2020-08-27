devtools::load_all("../../../dce")

library(tidyverse)


# parameters
target_dir <- strsplit(dirname(snakemake@output$graph_files[[1]]), "/")[[1]][[1]]
dir.create(file.path(target_dir, "csv_files"), recursive = TRUE)

database_list <- c("kegg", "biocarta", "pharmgkb")

# pathway meta information
df_info <- dce::get_pathway_info(database_list = database_list) %>%
  write_csv(file.path(target_dir, "pathway_info.csv"))

# make sure we retrieve (at least) all needed pathways
stopifnot(all(snakemake@params$pathways %in% str_remove(df_info$id, ":")))

# download pathways
options(Ncpus = 4)
dce::get_pathways(database_list = database_list) %>%
  purrr::map(function(x) {
    id_clean <- str_remove(x$id, ":")
    fname_csv <- file.path(target_dir, "csv_files", glue::glue("{id_clean}.csv")) # {x$database}_

    dce::graph2df(x$graph) %>%
      write_csv(fname_csv)
  }) %>%
  invisible
