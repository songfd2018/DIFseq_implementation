###########################################
# Differential abundance analysis by Milo #
###########################################
library(miloR)
library(SingleCellExperiment)
library(scater) # runPCA, runUMAP
library(dplyr) # distinct
library(ggplot2)
library(mclust) # adjustedRandIndex
library(VennDiagram)
library(RColorBrewer)

rm(list=ls())
proj <- "pancreas"
ver <- 1
method <- "Milo"

#############
# Load data #
#############
# Load data
y_obs <- read.table(paste0("RawCountData/count_data_",proj,"_v",ver,".txt"))
y_obs <- as.matrix(y_obs)

# Load dimension
dim <- read.table(paste0("RawCountData/dim_",proj,"_v",ver,".txt"))
dim <- unlist(dim)
N <- dim[1]
S <- dim[2]
G <- dim[3]
B <- dim[4]
Num_Treatment <- dim[5]
P <- dim[6]
BT_pair <- matrix(dim[6 + 1:(3*P)], byrow = TRUE, nrow = P)
colnames(BT_pair) <- c("Batch", "Treatment", "n_bt")

# Load metadata
metadata <- read.table(paste0("RawCountData/metadata_",proj,"_v",ver,".txt"), header = TRUE)
# colnames(metadata) <- c("Batch", "Treatment", "CellType")

# Load the gene list
gene_list <- unlist(read.table(paste0("RawCountData/gene_list_",proj,"_v",ver,".txt"),stringsAsFactors = F))

# Add row and column names for the raw count data matrix
rownames(y_obs) <- gene_list
colnames(y_obs) <- metadata$Sample

# indices of cells
btp_ind <- read.table(paste0("RawCountData/tbinfor_",proj,"_v",ver,".txt"))

# Build a SingleCellExperiment object
study_names <- unique(metadata$Study)

y_sub <- y_obs[, metadata$Study %in% study_names[3:4]]
metadata_sub <- metadata[metadata$Study %in% study_names[3:4], ]

sce_pancreas <- SingleCellExperiment(assays = list(counts = y_sub, logcounts = log1p(y_sub)),
                                     colData = metadata_sub)

# Load PCA and UMAP for visualization
PCA_embed <- read.csv(paste0("RawCountData/PCA_MNN_",proj,"_v",ver,".csv"))
UMAP_embed <- read.csv(paste0("RawCountData/UMAP_MNN_",proj,"_v",ver,".csv"))

rownames(PCA_embed) <- colnames(y_sub)
rownames(UMAP_embed) <- colnames(y_sub)

PCA_embed <- PCA_embed[,-1]
UMAP_embed <- UMAP_embed[,-1]

reducedDim(sce_pancreas, "pca") <- as.matrix(PCA_embed)
reducedDim(sce_pancreas, "umap") <- as.matrix(UMAP_embed)

############################
# Clustering for all cells #
############################
load("Output/pancreas_low_dimension_MNN_v1.RData")
# Create a Milo object
pan_milo <- Milo(fastCorrect)
pan_milo

# Construct KNN graph
pan_milo <- buildGraph(pan_milo, k = 30, d = 30, reduced.dim = "pca.corrected")
pan_milo <- makeNhoods(pan_milo, prop = 0.1, k = 30, d=30, refined = TRUE, reduced_dims = "pca.corrected")

# Counting cells in neighbourhoods
pan_milo <- countCells(pan_milo, meta.data = metadata, samples="Pair")
head(nhoodCounts(pan_milo))

## Convert batch info from integer to factor
pan_design <- metadata[,c("Pair", "Disease", "Study")]
pan_design <- distinct(pan_design)
rownames(pan_design) <- pan_design$Pair

# 3.6 computing neighbourhood connectivity
pan_milo <- calcNhoodDistance(pan_milo, d=30, reduced.dim = "pca.corrected")

# 3.7 testing
da_results <- testNhoods(pan_milo, design = ~ Disease, design.df = pan_design, reduced.dim = "UMAP")

# 5.1 Automatic grouping of neighbourhoods
pan_milo <- buildNhoodGraph(pan_milo)

## Find groups
da_results <- groupNhoods(pan_milo, da_results, da.fdr = 0.6, max.lfc.delta = 10)

# Assign cells to each Nhood Group
matr_cell_nhood <- nhoods(pan_milo)
cell_index <- summary(matr_cell_nhood)$i
Nhood_index <- summary(matr_cell_nhood)$j
cluster_milo <- rep(NA, nrow(metadata))
for(i in 1:N){
  nhoods_i <- which(cell_index==i)
  if(length(nhoods_i)>0){
    group_temp <- da_results$NhoodGroup[Nhood_index[nhoods_i]]
    freq <- table(group_temp)
    cluster_milo[i] <- names(freq[which.max(freq)])
  }
}

ARI_MiloR <- adjustedRandIndex(cluster_milo,metadata$CellType)
write.table(ARI_MiloR, file = paste0("Output/ARI_",method,"_",proj,"_v",ver,".txt"), quote = FALSE)

##########################################################
# Differential abundance testing on the last two batches #
##########################################################
# Refer to https://rawcdn.githack.com/MarioniLab/miloR/7c7f906b94a73e62e36e095ddb3e3567b414144e/vignettes/milo_gastrulation.html#1_Load_data
# Create a Milo object
pancreas_milo <- Milo(sce_pancreas)
pancreas_milo

# Construct KNN graph
pancreas_milo <- buildGraph(pancreas_milo, k = 30, d = 30, reduced.dim = "pca")
pancreas_milo <- makeNhoods(pancreas_milo, prop = 0.1, k = 30, d=30, refined = TRUE, reduced_dims = "pca")
plotNhoodSizeHist(pancreas_milo)

# Counting cells in neighbourhoods
metadata_sub$Pair <- paste0(metadata_sub$Study, "_", metadata_sub$Disease)

pancreas_milo <- countCells(pancreas_milo, meta.data = metadata_sub, samples="Pair")
# head(nhoodCounts(pancreas_milo))
# dim(nhoodCounts(pancreas_milo))

# Defining experimental design
pancreas_design <- metadata_sub[,c("Pair", "Disease", "Study")]

## Convert batch info from integer to factor
pancreas_design$Study <- factor(pancreas_design$Study, levels = study_names[3:4]) 
pancreas_design$Disease <- as.factor(pancreas_design$Disease) 
pancreas_design <- distinct(pancreas_design)
rownames(pancreas_design) <- pancreas_design$Pair

# 3.6 computing neighbourhood connectivity
pancreas_milo <- calcNhoodDistance(pancreas_milo, d=30, reduced.dim = "pca")

# 3.7 testing
da_results <- testNhoods(pancreas_milo, design = ~ Disease, design.df = pancreas_design, reduced.dim = "umap")
# head(da_results)
# 
# da_results %>%
#   arrange(SpatialFDR) %>%
#   head() 

# take the threshold as sptial FDR < 0.1
sum(da_results$SpatialFDR < 0.1)

##############################
# Draw UMAP of neighborhoods #
##############################
# Histogram of p values
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)

# Volcano plot
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) +
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

## Plot neighbourhood graph
pancreas_milo <- buildNhoodGraph(pancreas_milo)
nh_graph_pl <- plotNhoodGraphDA(pancreas_milo, da_results, layout="umap", alpha=0.1) 

type <- "UMAP_neighborhoods"
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,".pdf"),width = 12, height = 8)
nh_graph_pl
dev.off()

#######################################################
# Draw UMAP of each individual cell via average logFC #
#######################################################
# Compute average logFC for each cell
matr_cell_nhood <- nhoods(pancreas_milo)
cell_index <- summary(matr_cell_nhood)$i
Nhood_index <- summary(matr_cell_nhood)$j

N <- nrow(metadata_sub)
average_logFC <- rep(0, N)
for(i in 1:N){
  nhoods_i <- which(cell_index==i)
  if(length(nhoods_i)>0){
    
    average_logFC[i] <- mean(da_results$logFC[Nhood_index[nhoods_i]])
    # group_temp <- da_results$NhoodGroup[Nhood_index[nhoods_i]]
    # freq <- table(group_temp)
    # cluster_milo[i] <- names(freq[which.max(freq)])
  }
}

metadata_sub$AvelogFC <- average_logFC
metadata_sub$UMAP1 <- UMAP_embed[,1]
metadata_sub$UMAP2 <- UMAP_embed[,2]

method <- "Milo"
type <- "UMAP"
feature1 <- "UMAP1"
feature2 <- "UMAP2"
color_by <- "AvelogFC"
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,".pdf"),width = 12, height = 8)
p <- ggplot(metadata_sub, aes_string(x = feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_gradient2(low = "#FF80FF", mid = "#e0e0e0",high = "#80FFFF", midpoint = 0) +
  xlab(feature1) + ylab(feature2) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333", 
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue", 
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2)) +
  labs(color = "Ave. LFC")
print(p)
dev.off()

