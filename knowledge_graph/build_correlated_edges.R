library(dplyr)
setwd("~/Desktop/Research/")

edges <- read.csv("./unweightedControl/edges.csv")
nodes <- read.csv("./unweightedControl/nodes.csv")
kg <- read.csv("./unweightedControl/kg.csv")
ppi <- edges %>% filter(display_relation == "ppi")

cor_matrix <- readRDS("./KGAugmentation/phyloP_correlation.rds")
# cor_matrix_hs <- readRDS("KGAugmentation/phyloP_correlation_HS.rds")
active_cormatrix <- cor_matrix

ppi$relation <- "phyloP"
ppi$display_relation <- "phyloP"
ppi$x_name <- nodes$node_name[match(ppi$x_index, nodes$node_index)]
ppi$y_name <- nodes$node_name[match(ppi$y_index, nodes$node_index)]

ppi_intersection <- ppi %>% filter(x_name %in% intersect(x_name, colnames(active_cormatrix))
                                   & y_name %in% intersect(y_name, colnames(active_cormatrix)))

# create phyloP correlation edges alongside existing ppi edges
ppi_intersection$value <- mapply(function(row, col) active_cormatrix[row, col],
                 ppi_intersection$x_name, ppi_intersection$y_name)
ppi_intersection <- ppi_intersection %>% filter(!is.na(value))

ppi <- ppi %>% select(-x_name, -y_name)
ppi_intersection <- ppi_intersection %>% select(-x_name, -y_name)

correlated.edges <- ppi_intersection
correlated.kg <- data.frame(
  relation="phyloP",
  display_relation="phyloP",
  x_index=nodes$node_index[match(correlated.edges$x_index, nodes$node_index)],
  x_id=nodes$node_id[match(correlated.edges$x_index, nodes$node_index)],
  x_type=nodes$node_type[match(correlated.edges$x_index, nodes$node_index)],
  x_name=nodes$node_name[match(correlated.edges$x_index, nodes$node_index)],
  x_source=nodes$node_source[match(correlated.edges$x_index, nodes$node_index)],
  y_index=nodes$node_index[match(correlated.edges$x_index, nodes$node_index)],
  y_id=nodes$node_id[match(correlated.edges$y_index, nodes$node_index)],
  y_type=nodes$node_type[match(correlated.edges$y_index, nodes$node_index)],
  y_name=nodes$node_name[match(correlated.edges$y_index, nodes$node_index)],
  y_source=nodes$node_source[match(correlated.edges$y_index, nodes$node_index)],
  value=correlated.edges$value
)
