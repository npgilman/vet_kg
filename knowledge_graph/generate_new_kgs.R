library(dplyr)
setwd("~/Desktop/Research/")

edges <- read.csv("./unweightedControl/edges.csv")
nodes <- read.csv("./unweightedControl/nodes.csv")
kg <- read.csv("./unweightedControl/kg.csv")

### add weights of 1 to the edges and kg
edges$value <- 1
kg$value <- 1
write.csv(edges, file=paste0("./weightedControl/kg.csv"), row.names=FALSE)
write.csv(kg, file=paste0("./weightedControl/edges.csv"), row.names=FALSE)

### for this section, must have the following variables generated from build_self_loops.R
###### self_loops.edges
###### self_loops.kg
output.sl.edges <- rbind(self_loops.edges, edges)
output.sl.kg <- rbind(self_loops.kg, kg)
write.csv(output.sl.kg, file=paste0("./weightedSelfLoops/kg.csv"), row.names=FALSE)
write.csv(output.sl.edges, file=paste0("./weightedSelfLoops/edges.csv"), row.names=FALSE)
View(output.sl.kg)
nrow(output.sl.kg) - nrow(kg)

### for this section, must have the following variables generated from build_correlated_edges.R
###### correlated.edges
###### correlated.kg
output.c.edges <- rbind(correlated.edges, edges)
output.c.kg <- rbind(correlated.kg, kg)
write.csv(output.c.kg, file=paste0("./weightedCorrelation/kg.csv"), row.names=FALSE)
write.csv(output.c.edges, file=paste0("./weightedCorrelation/edges.csv"), row.names=FALSE)
View(output.c.kg)
nrow(output.c.kg) - nrow(kg)

kg<-read.csv("./weightedCorrelation/kg.csv")
View(kg)
