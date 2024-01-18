library(batchelor)
library(scater) # runPCA function

rm(list=ls())
proj <- "covid"
ver <- 1

#############################################
# Load read count and dimension information #
#############################################
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

########################
#### MNN correction ####
########################
method <- "MNN"

# Refer to https://bioconductor.org/packages/release/bioc/vignettes/batchelor/inst/doc/correction.html
# Create the separate singlecellExperiments object for each batch
sce <- list()
for(b in 1:B){
  
  cell_index <- metadata$Batch == b
  sce[[b]] <- SingleCellExperiment(assays = list(counts = y_obs[,cell_index],
                                                 logcounts = log1p(y_obs[,cell_index])),
                                   colData = metadata[cell_index,])
}

# Combine data sets from 6 batches together without correction
combined <- correctExperiments(Batch1 = sce[[1]], Batch2 = sce[[2]],
                               PARAM=NoCorrectParam())

# batch effect correction by fastMNN()
fastCorrect <- fastMNN(combined, batch = combined$batch)
fastCorrect <- runPCA(fastCorrect, dimred = "corrected", name = "pca.corrected")

# Generate UMAP
fastCorrect <- runUMAP(fastCorrect, dimred = "pca.corrected", pca = 30,
                       name = "UMAP", n_neighbors = 100)

PCA_MNN <- reducedDim(fastCorrect, type = "pca.corrected")
UMAP_MNN <- reducedDim(fastCorrect, type = "UMAP")

n.celltype <- length(unique(metadata$CellType))

#############################################
# Store the dimension reduction information #
#############################################
PCA_output <- PCA_MNN
UMAP_output <- UMAP_MNN
metadata_output <- metadata

# metadata_output$Cluster <- as.numeric(factor(metadata_output$DIFseq_cluster))
metadata_output$Cluster_l1 <- as.numeric(factor(metadata_output$Seurat_Level1))
metadata_output$Cluster_l2 <- as.numeric(factor(metadata_output$Seurat_Level2))
metadata_output$Batch <- paste0("Batch",metadata_output$Batch)

# rownames(PCA_output) <- paste0("Sample",metadata$sample, "_Cell", rep(1:ns[1], S))
colnames(PCA_output) <- paste0("PC",1:ncol(PCA_MNN))

# rownames(UMAP_output) <- rownames(PCA_output)
colnames(UMAP_output) <- c("UMAP1", "UMAP2")

write.csv(PCA_output, file = paste0("RawCountData/PCA_",method,"_",proj,"_v",ver,".csv"), quote = FALSE)
write.csv(UMAP_output, file = paste0("RawCountData/UMAP_",method,"_",proj,"_v",ver,".csv"), quote = FALSE)
write.csv(metadata_output, file = paste0("RawCountData/metadata_",proj,"_v",ver,".csv"), quote = FALSE)

save.image(paste0("Output/",method,"_analysis_",proj,"_v",ver,".RData"))