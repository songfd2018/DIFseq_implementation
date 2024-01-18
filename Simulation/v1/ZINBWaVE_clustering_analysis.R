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
ver <- 1
method <- "ZINBWaVE"

#############
# Load data #
#############

setwd(paste0("Simulation/v",ver))
load(paste0("simulation_countdata_v",ver,".RData"))

# setwd("/home/songfangda/DIFseq/Model_Comparison/ZINBWaVE/")
# load(paste0("../simulation_countdata_v",ver,".RData"))

sce <- SingleCellExperiment(assays = list(counts = y))

# Build three dummy variables to correct for batch effects 
sce$Batch2 <- metadata$batch == 2
sce$Batch3 <- metadata$batch == 3
sce$Batch4 <- metadata$batch == 4
sce$Cond <- metadata$condition

###################
# Apply ZINB-WaVE #
###################
ZINBWaVE_start <- Sys.time()
simulation_cov <- zinbwave(sce, K=2, X="~Batch2 + Batch3 + Batch4", 
                           epsilon=G, BPPARAM=BiocParallel::SerialParam())
ZINBWaVE_end <- Sys.time()

print(paste0("It takes ", difftime(ZINBWaVE_end, ZINBWaVE_start, units = "mins"), " to run ZINBWaVE for simulation_v",ver,"."))

# save.image(paste0("Output/",method,"_analysis_",proj,"_v",ver,".RData"))

########################
# Clustering by Seurat #
########################
# Convert the SingleCellExperiment object to a Seurat object
seu <- as.Seurat(x = simulation_cov, counts = "counts", data = "counts")

# Finally, we can use the zinbwave factors for cluster analysis.
seu <- FindNeighbors(seu, reduction = "zinbwave",
                     dims = 1:2)
seu <- FindClusters(object = seu, resolution = 0.2)

ARI <- adjustedRandIndex(seu$seurat_clusters, metadata$celltype)
write.table(ARI, file = paste0("Output/ARI_",method,"_",proj,"_v",ver,".txt"), quote = FALSE)
