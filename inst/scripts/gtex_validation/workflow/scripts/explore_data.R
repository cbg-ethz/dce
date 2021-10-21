###
# Plot pathway under consideration.
###


library(tidyverse)
library(cowplot)
devtools::load_all("../../..")

View(dce::df_pathway_statistics %>%
       arrange(desc(node_num)))
pathways <- get_pathways(pathway_list = list(kegg = c("Pathways in cancer")))
pathway <- pathways[[1]]$graph
plot_network(as(pathway, 'matrix'), visualize_edge_weights = FALSE)
