library(dplyr)
library(stringr)
library(tidyselect)
setwd("~/Desktop/Research/")

### loading primekg unweighted control for EDA
{
  nodes <- read.csv("unweightedControl/nodes.csv")
  drugs <- nodes %>% filter(node_type == "drug")
  proteins <- nodes %>% filter(node_type == "gene/protein")
  View(proteins)
  
  # load edges/kg and filter for drug_protein interactions
  edges <- read.csv("unweightedControl/edges.csv")
  drug_protein_edges <- edges %>% filter(str_detect(relation, "drug_protein"))
  kg <- read.csv("unweightedControl/kg.csv")
  drug_protein_kg <- kg %>% filter(str_detect(relation, "drug_protein"))
}

### exploring primekg drug features from drug_bank and drug_central
{
drug_features <- read.csv("unweightedControl/drug_features.csv")
drugwnames <- drug_features %>% mutate(name = nodes$node_name[node_index + 1]) %>% relocate(name)
View(drugwnames)
}

### loading drug central interaction tsv
{
drugcentral <- read.csv("drug_central/drug.target.interaction.tsv", sep="\t")
dgo <- drugcentral %>% select(DRUG_NAME, GENE, ORGANISM) %>% mutate(GENE = toupper(GENE))
# # num unique genes in drugcentral tsv
# length(unique(tolower(drugcentral$GENE)))
# # num unique drugs in drugcentral tsv
# length(unique(tolower(drugcentral$DRUG_NAME)))
}

# comparing the intersection of:
#      drug central :  2140 unique genes | 2587 unique drugs
#      primeKG      : 27671 unique genes | 7957 unique drugs
#      intersection :  1623 unique genes | 2126 unqiue drugs
#                   
## note: primeKG includes genes and proteins as a single node type
## note: only a 5% increase in drugs added to primeKG
{
#drugs
drug_intersection <- intersect(unique(tolower(drugcentral$DRUG_NAME)), unique(tolower(drugs$node_name)))
length(unique(tolower(drugs$node_name))) # unique primeKG
length(unique(tolower(drugcentral$DRUG_NAME))) # unique drug central
length(drug_intersection)

#genes
gene_intersection <- intersect(unique(tolower(drugcentral$GENE)), unique(tolower(proteins$node_name)))
length(unique(tolower(proteins$node_name)))
length(unique(tolower(drugcentral$GENE)))
length(gene_intersection)
}

# subset drugcentral data for later creation of new drug nodes and edges
{
new_drugs <- drugcentral %>% filter(! (tolower(DRUG_NAME) %in% drug_intersection) )
new_drugs_old_genes_intersection <- intersect(unique(tolower(nodes$node_name)), unique(tolower(new_drugs$GENE)))
# 641 / 769 unique genes associated with new_drugs are already in primeKG. 
# omit drugs associated with the missing 127 genes for now - this removes 54 drugs from being added
new_drugs_old_genes <- new_drugs %>% filter( tolower(GENE) %in% new_drugs_old_genes_intersection)
new_drugs_new_genes <- new_drugs %>% filter(! (tolower(GENE) %in% new_drugs_old_genes_intersection) )
# pie(table(new_drugs_old_genes$ORGANISM))
# pie(table(new_drugs_new_genes$ORGANISM))
# sort(table(new_drugs_new_genes$ORGANISM)) # see a table of what organisms contain newly introduced genes 

# explore how much new of the genes that come from new drugs from non-majority species are retained
{
non_majority <- new_drugs %>% filter(! ((ORGANISM == "Homo sapiens") | 
                                     (ORGANISM == "Rattus norvegicus") | 
                                     (ORGANISM == "Mus musculus") |
                                     (ORGANISM == "Bos taurus")))

length(unique(tolower(non_majority$GENE)))
length(intersect(unique(tolower(nodes$node_name)), unique(tolower(non_majority$GENE))))
}



# build new drug nodes in format of primeKG
{
new_drugs.nodes <- data.frame(matrix(ncol=5, nrow=0))
colnames(new_drugs.nodes) <- c("node_index", "node_id", "node_type", "node_name", "node_source")

node_index <- 129374
node_type <- "drug"
node_source <- "DrugCentral"

unique_drugs <- unique(tolower(new_drugs_old_genes$DRUG_NAME))
length(unique_drugs)
for (i in seq(1,length(unique_drugs))) {
  new_drug.row <- data.frame(node_index = node_index + 1,
                             node_id = sprintf("DC%05d", i),
                             node_type = node_type,
                             node_name = unique_drugs[i],
                             node_source = node_source)

  new_drugs.nodes <- rbind(new_drugs.nodes, new_drug.row)
  node_index <- node_index + 1
}
}

# build new drugs edges in the format of primeKG
{
new_drugs.edges <- data.frame(matrix(ncol=4, nrow=0))
colnames(new_drugs.edges) <- c("relation", "display_relation", "x_index", "y_index")

prev_drug <- ""
index <- 0
gene_index_match <- nodes$node_index[match(tolower(new_drugs_old_genes$GENE), tolower(nodes$node_name))]

# make sure to run INDEX <- 0 before running below loop
for (i in seq(1,nrow(new_drugs_old_genes))) {
  if (prev_drug != new_drugs_old_genes$DRUG_NAME[i]) {
    print(paste0("old: ", prev_drug, " current: ", new_drugs_old_genes$DRUG_NAME[i]))
    index <- index + 1
    prev_drug <- new_drugs_old_genes$DRUG_NAME[i]
  }
  
  # skip duplicate edges due to multiple organisms with same genes
  # only done by checking previous edges
  new_drugs.edges$x_index[i - 1] == gene_index_match[i]
  new_drugs.edges
  if (nrow(new_drugs.edges) > 1 && (new_drugs.edges$x_index[i - 1] == gene_index_match[i])) {next}
  print(new_drugs.nodes$node_index[index])
  new_drugs.row1 <- data.frame(relation= "drug_protein",
                               display_relation= "N/A",
                               x_index= new_drugs.nodes$node_index[index],
                               y_index= nodes$node_index[gene_index_match[i]])
  # new_drugs.row2 <- data.frame(relation= "drug_protein",
  #                             display_relation= "N/A",
  #                             x_index= nodes$node_index[gene_index_match[i]],
  #                             y_index= new_drugs.nodes$node_index[index])
  
  new_drugs.edges <- rbind(new_drugs.edges, new_drugs.row1)
}

# edges <- rbind(new_drugs.edges, edges)
}

# 135 VET - 216 HUMVET Drugs
{
# for the below line to work, nodes must be rbind with new_drugs.nodes
nodes <- rbind(new_drugs.nodes, nodes)
new_drugs.kg <- data.frame(relation="drug_protein",
                            display_relation="N/A",
                            x_index=nodes$node_index[match(new_drugs.edges$x_index, nodes$node_index)],
                            x_id=nodes$node_id[match(new_drugs.edges$x_index, nodes$node_index)],
                            x_type=nodes$node_type[match(new_drugs.edges$x_index, nodes$node_index)],
                            x_name=nodes$node_name[match(new_drugs.edges$x_index, nodes$node_index)],
                            x_source=nodes$node_source[match(new_drugs.edges$x_index, nodes$node_index)],
                            y_index=nodes$node_index[match(new_drugs.edges$y_index, nodes$node_index)],
                            y_id=nodes$node_id[match(new_drugs.edges$y_index, nodes$node_index)],
                            y_type=nodes$node_type[match(new_drugs.edges$y_index, nodes$node_index)],
                            y_name=nodes$node_name[match(new_drugs.edges$y_index, nodes$node_index)],
                            y_source=nodes$node_source[match(new_drugs.edges$y_index, nodes$node_index)])

# View(new_drugs.kg)
}

sum(is.na(rbind(new_drugs.kg, kg)))
write.csv(rbind(new_drugs.edges, edges), file=paste0("./drug_central/vet_drugs/edges.csv"), row.names=FALSE)
write.csv(rbind(new_drugs.kg, kg), file=paste0("./drug_central/vet_drugs/kg.csv"), row.names=FALSE)
write.table(rbind(new_drugs.nodes, nodes), file=paste0("./drug_central/vet_drugs/node.csv"), sep="\t", row.names=FALSE)

# viewing data format for bidirectional edges with example [Copper <--> F8]
View(
  kg %>% 
    filter(x_name == "Copper" | y_name == "Copper") %>%
    filter(x_name == "F8" | y_name == "F8") %>%
    filter(relation == "drug_protein")
)
