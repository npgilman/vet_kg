library(dplyr)
setwd("~/Desktop/Research")

edges <- read.csv("primeKG/edges.csv")
nodes <- read.csv("primeKG/nodes.csv") %>% filter(node_type == "gene/protein")
kg <- read.csv("primeKG/kg.csv")

phyloP <- read.csv("features.tsv", sep="\t") %>%
  filter(GeneID %in% nodes$node_name) %>% 
  filter(!is.na(mean_Homo_sapiens))

matched_indices <- match(phyloP$GeneID, nodes$node_name)
self_loops.nodes <- data.frame(node_index=matched_indices - 1,
                               node_id=nodes$node_id[matched_indices],
                               node_type="gene/protein",
                               node_name=phyloP$GeneID,
                               node_source="human gtf file")

self_loops.edges <- data.frame(relation="phyloP",
                               display_relation="phyloP",
                               x_index=self_loops.nodes$node_index,
                               y_index=self_loops.nodes$node_index,
                               value=phyloP$mean_Homo_sapiens)

self_loops.kg <- data.frame(relation="phyloP",
                            display_relation="phyloP",
                            x_index=self_loops.nodes$node_index,
                            x_id=self_loops.nodes$node_id,
                            x_type=self_loops.nodes$node_type,
                            x_name=self_loops.nodes$node_name,
                            x_source=self_loops.nodes$node_source,
                            y_index=self_loops.nodes$node_index,
                            y_id=self_loops.nodes$node_id,
                            y_type=self_loops.nodes$node_type,
                            y_name=self_loops.nodes$node_name,
                            y_source=self_loops.nodes$node_source,
                            value=phyloP$mean_Homo_sapiens)

