#' Dataframe containing meta-information of pathways in database
#'
#' @param query_species For which species
#' @param database_list Which databases to query. Query all if `NULL`.
#' @param include_network_statistics Compute some useful statistics
#'        per pathway. Takes longer!
#' @import graphite graph glue purrr logger
#' @importFrom dplyr pull
#' @importFrom rlang .data
#' @export
#' @return data frame with pathway meta information
#' @examples
#' head(get_pathway_info(database_list = c("kegg")))
get_pathway_info <- function(
    query_species = "hsapiens", database_list = NULL,
    include_network_statistics = FALSE
) {
    if (is.null(database_list)) {
        database_list <- graphite::pathwayDatabases() %>%
            dplyr::filter(.data$species == query_species) %>%
            pull(.data$database)
    }

    database_list %>%
        purrr::map_dfr(function(database) {
            logger::log_info("Processing {database}")
            db <- graphite::pathways(query_species, database)

            purrr::map_dfr(as.list(db), function(pw) {
                tmp <- data.frame(
                    database = database,
                    pathway_id = graphite::pathwayId(pw),
                    pathway_name = graphite::pathwayTitle(pw)
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
#' @param pathway_list List mapping database name to
#' vector of pathway names to download
#' @import graphite glue purrr org.Hs.eg.db logger
#' @importFrom rlang .data
#' @export
#' @return list of pathways
#' @examples
#' pathways <- get_pathways(
#'   pathway_list = list(kegg = c(
#'     "Protein processing in endoplasmic reticulum"
#'   ))
#' )
#' plot_network(as(pathways[[1]]$graph, "matrix"))
get_pathways <- function(
    query_species = "hsapiens", database_list = NULL,
    remove_empty_pathways = TRUE,
    pathway_list = NULL
) {
    if (is.null(database_list)) {
        # we need to figure out which databases to use
        if (is.null(pathway_list)) {
            # user has no preference, use all
            database_list <- graphite::pathwayDatabases() %>%
                filter(.data$species == query_species) %>%
                pull(.data$database)
        } else {
            # user wants certain pathways, only use respective databases
            database_list <- names(pathway_list)
        }
    }

    database_list %>%
        purrr::map(function(database) {
            logger::log_info("Processing {database}")

            db <- graphite::pathways(query_species, database)

            if (!is.null(pathway_list) && !is.null(pathway_list[[database]])) {
                db <- db[pathway_list[[database]]]
            }

            db_symbol <- graphite::convertIdentifiers(db, "SYMBOL")

            purrr::map(as.list(db_symbol), function(pw) {
                # remove "SYMBOL:" prefix
                graph <- graphite::pathwayGraph(pw, which = "proteins")

                if (length(nodes(graph)) != 0) {
                    nodes(graph) <- vapply(nodes(graph), function(x) {
                        strsplit(x, ":")[[1]][[2]]
                    }, FUN.VALUE = character(1), USE.NAMES = FALSE)
                }

                if (remove_empty_pathways & length(nodes(graph)) == 0) {
                    return(NULL)
                }

                list(
                    database = database,
                    pathway_id = graphite::pathwayId(pw),
                    pathway_name = graphite::pathwayTitle(pw),
                    graph = graph
                )
            }) %>%
                unname
        }) %>%
        unlist(recursive = FALSE, use.names = FALSE) %>%
        purrr::compact()
}
