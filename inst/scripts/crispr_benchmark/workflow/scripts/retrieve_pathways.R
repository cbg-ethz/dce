###
# Retrieve pathway data from all suppored databases.
###


library(tidyverse)

devtools::load_all("../../../")


# parameters
target_dir <- snakemake@output$pathway_dir
dir.create(file.path(target_dir, "csv_files"), recursive = TRUE)

database_list <- snakemake@config$databases

# pathway meta information
df_info <- dce::get_pathway_info(database_list = database_list) %>%
  write_csv(file.path(target_dir, "pathway_info.csv"))

# download pathways
options(Ncpus = 4)
dce::get_pathways(database_list = database_list) %>%
  purrr::map(function(x) {
    id_clean <- str_remove(x$pathway_id, ":")
    fname_csv <- file.path(target_dir, "csv_files", glue::glue("{id_clean}.csv")) # {x$database}_

    dce::graph2df(x$graph) %>%
      write_csv(fname_csv)
  }) %>%
  invisible
