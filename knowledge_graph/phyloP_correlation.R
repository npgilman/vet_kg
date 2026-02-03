library(dplyr)
library(tibble)
library(WGCNA)
library(parallel)
library(BiocGenerics)
library(pheatmap)
detectCores()
setwd("~/Desktop/Research/KGAugmentation/")

### TxGNN Data
nodes <- read.csv("../unweightedControl/nodes.csv")
genes <- filter(nodes, nodes$node_type == "gene/protein")

### PhyloP Data
gene_feat <- read.csv("../features.tsv", sep="\t")
gene_feat_HS <- gene_feat %>%
  select("GeneID" | contains("Homo")) %>% 
  filter(! if_any(everything(), is.na))

GeneID <- gene_feat$GeneID
gene_feat$GeneID <- NULL
GeneID_HS <- gene_feat_HS$GeneID
gene_feat_HS$GeneID <- NULL

t_gf <- t(gene_feat)
t_gfhs <- t(gene_feat_HS)
colnames(t_gf) <- GeneID
colnames(t_gfhs) <- GeneID_HS
# View(t_gfhs)

### Sort IQR for homosapiens and grab the a subset of variable genes
IQR_hs <- gene_feat_HS$IQR_Homo_sapiens
names(IQR_hs) <- GeneID_HS
IQR_hs <- sort(IQR_hs, decreasing=TRUE)
top1000 <- IQR_hs[1:1000]
bottom1000 <- tail(IQR_hs, 1000)
gene_subset <- c(unlist(names(top1000)), unlist(names(bottom1000)))
gene_subset <- names(IQR_hs[1:1000])
# gene_subset
# head(top1000)

iqr1000 <- IQR_hs[seq(1, length(IQR_hs), by=18)]
gene_subset <- names(iqr1000)
# head(gene_subset)

### Calculate cor matrix
enableWGCNAThreads(nThreads=14)
# cor_matrix <- cor(t_gf, use = "pairwise.complete.obs", nThreads=14)
# saveRDS(cor_matrix, file="phyloP_correlation.rds")

cor_matrix <- readRDS("phyloP_correlation.rds")
# head(cor_matrix)


# cor_matrix_hs <- cor(t_gfhs, nThreads=14)
# saveRDS(cor_matrix_hs, file="phyloP_correlation_HS.rds")
cor_matrix_hs <- readRDS("phyloP_correlation_HS.rds")
# head(cor_matrix_hs)

# cor_vals <- as.list(cor_matrix_hs)
# summary(cor_vals)

### Display results
# View(cor_matrix)
subset_mat <- cor_matrix[gene_subset, gene_subset]
subset_mat_hs <- cor_matrix_hs[gene_subset, gene_subset]

# heatmap(subset_mat, main="corr matrix (top1000 + bottom1000)", col =
# colorRampPalette(c("blue", "white", "red"))(100))

heatmap(subset_mat_hs,
        main="corr matrix homosapiens",
        labX=NULL,
        labY=NULL,
        col = colorRampPalette(c("blue", "white", "red"))(100))

### Extract correlation as a edge
length(GeneID)
gg <- intersect(genes$node_name, GeneID)


vals <- c()
count <- 0
for (i in 1:length(gg)) {
  for (j in i:length(gg)) {
    count <- count + 1
    vals <- c(vals, cor_matrix[i, j])
  }
#  print(count)
}

