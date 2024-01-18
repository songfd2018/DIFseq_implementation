##################################### 
# Cell type clustering by ZINB-WaVE #
#####################################
library(zinbwave)
library(SingleCellExperiment)

library(mclust)
library(Seurat)
library(viridis)
library(ggplot2)
# library(mclust)

rm(list=ls())
proj <- "simulation"
ver <- 3
method <- "ZINBWaVE"

#############
# Load data #
#############

setwd(paste0("Simulation/v",ver))
y <- read.table(file = paste0("RawCountData/count_data_",proj,"_v",ver,".txt"))
y <- as.matrix(y)

dim <- read.table(paste0("RawCountData/dim_",proj,"_v",ver,".txt"))
dim <- unlist(dim)
N <- dim[1]
S <- dim[2]
G <- dim[3]
B <- dim[4]
Num_Treatment <- dim[5]
P <- dim[6]

BT_pair <- matrix(dim[5 + 1:(3*P)], byrow = TRUE, nrow = P) 
colnames(BT_pair) <- c("Batch", "Condition", "n_bt")

# Load metadata
metadata <- read.table(paste0("RawCountData/metadata_",proj,"_v",ver,".txt"))

# load(paste0("../simulation_countdata_v",ver,".RData"))

sce <- SingleCellExperiment(assays = list(counts = y))

# Build three dummy variables to correct for batch effects 
sce$Batch2 <- metadata$batch == 2
sce$Batch3 <- metadata$batch == 3
sce$Batch4 <- metadata$batch == 4
sce$Batch5 <- metadata$batch == 5
sce$Batch6 <- metadata$batch == 6 
sce$Cond <- metadata$condition

##############################
# Joint matrix factorization #
##############################
ZINBWaVE_start <- Sys.time()
simulation_cov <- zinbwave(sce, K=2, X="~Batch2 + Batch3 + Batch4 + Batch5 + Batch6", 
                  epsilon=G, BPPARAM=BiocParallel::SerialParam())
ZINBWaVE_end <- Sys.time()
print(paste0("It takes ", difftime(ZINBWaVE_end, ZINBWaVE_start, units = "mins"), " to run ZINBWaVE for simulation_v",ver,"."))

###########################
# Do clustering by Seurat #
###########################
# Convert the SingleCellExperiment object to a Seurat object
rownames(simulation_cov) <- paste0("Genes_",1:G)
seu <- as.Seurat(x = simulation_cov, counts = "counts", data = "counts")

# Finally, we can use the zinbwave factors for cluster analysis.
seu <- FindNeighbors(seu, reduction = "zinbwave",
                     dims = 1:2)
seu <- FindClusters(object = seu, resolution = 0.2)

ARI <- adjustedRandIndex(seu$seurat_clusters, metadata$celltype)
write.table(ARI, file = paste0("Output/ARI_",method,"_",proj,"_v",ver,".txt"), quote = FALSE)
