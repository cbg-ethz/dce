#' Dataframe containing meta-information of pathways in database
#'
#' @param query_species For which species
#' @param database_list Which databases to query. Query all if `NULL`.
#' @param include_network_statistics Compute some useful statistics
#'        per pathway. Takes longer!
#' @import graphite graph glue purrr
#' @importFrom dplyr pull
#' @export
#' @return data frame with pathway meta information
#' @examples
#' print('example needed')
get_pathway_info <- function(
  query_species = "hsapiens", database_list = NULL,
  include_network_statistics = FALSE
) {
  if (is.null(database_list)) {
    database_list <- graphite::pathwayDatabases() %>%
      filter(species == query_species) %>%
      pull(database)
  }

  database_list %>%
    purrr::map_dfr(function(database) {
      print(glue::glue("Processing {database}"))
      db <- graphite::pathways(query_species, database)

      purrr::map_dfr(as.list(db), function(pw) {
        tmp <- data.frame(
          database = database,
          id = graphite::pathwayId(pw),
          name = graphite::pathwayTitle(pw)
        )

        if (include_network_statistics) {
          graph <- graphite::pathwayGraph(pw, which = "proteins")

          tmp$node_num <- graph::numNodes(graph)
          tmp$edge_num <- graph::numEdges(graph)
        }

        return(tmp)
      })
    })
}


#' Easy pathway network access
#'
#' @param query_species For which species
#' @param database_list Which databases to query. Query all if `NULL`.
#' @param remove_empty_pathways Discard pathways without nodes
#' @import graphite glue purrr
#' @export
#' @return list of pathways
#' @examples
#' print('example needed')
get_pathways <- function(
  query_species = "hsapiens", database_list = NULL,
  remove_empty_pathways = TRUE
) {
  if (is.null(database_list)) {
    database_list <- graphite::pathwayDatabases() %>%
      filter(species == query_species) %>%
      pull(database)
  }

  database_list %>%
    purrr::map(function(database) {
      print(glue::glue("Processing {database}"))

      db <- graphite::pathways(query_species, database)
      db_symbol <- graphite::convertIdentifiers(db, "SYMBOL")

      graph_list <- purrr::map(as.list(db_symbol), function(pw) {
        # remove "SYMBOL:" prefix
        graph <- graphite::pathwayGraph(pw, which = "proteins")
        if (length(nodes(graph)) != 0) {
          nodes(graph) <- vapply(nodes(graph), function(x) {
            strsplit(x, ":")[[1]][[2]]
          }, USE.NAMES = FALSE)
        }

        if (remove_empty_pathways & length(nodes(graph)) == 0) {
          return(NULL)
        }

        list(
          database = database,
          id = graphite::pathwayId(pw),
          name = graphite::pathwayTitle(pw),
          graph = graph
        )
      }) %>%
        unname
    }) %>%
    unlist(recursive = FALSE, use.names = FALSE) %>%
    purrr::compact()
}
