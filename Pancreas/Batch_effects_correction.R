library(Seurat)
library(ggplot2)
library(batchelor)
library(scater) # runPCA function
library(viridis)
library(mclust)
library(RColorBrewer)

#############################################
# Load read count and dimension information #
#############################################
rm(list=ls())
proj <- "pancreas"
ver <- 1

#####################
# Seurat correction #
#####################
method <- "Seurat"

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
study_names <- unique(metadata$Study)

# Load the gene list
gene_list <- unlist(read.table(paste0("RawCountData/gene_list_",proj,"_v",ver,".txt"),stringsAsFactors = F))

# Add row and column names for the raw count data matrix
rownames(y_obs) <- gene_list
colnames(y_obs) <- metadata$Sample

# indices of cells
btp_ind <- read.table(paste0("RawCountData/tbinfor_",proj,"_v",ver,".txt"))

pancreas.list <- NULL


colnames(y_obs) <- metadata$Sample
rownames(y_obs) <- gene_list
rownames(metadata) <- metadata$Sample

for(b in 1:B){
  
  cell_index <- metadata$Study == study_names[b]
  pancreas.list[[b]] <- CreateSeuratObject(y_obs[,cell_index], 
                                             project = study_names[b], 
                                             meta.data = metadata[cell_index,])
}

pancreas.list <- lapply(X = pancreas.list, FUN = function(x) {
  x <- NormalizeData(x)
})

# Here, we set the features of interest as the gene set analyzed in DIFseq
pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, anchor.features = gene_list)
pancreas.combined <- IntegrateData(anchorset = pancreas.anchors)
DefaultAssay(pancreas.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
pancreas.combined <- ScaleData(pancreas.combined, verbose = FALSE)
pancreas.combined <- RunPCA(pancreas.combined, npcs = 50, verbose = FALSE)

ElbowPlot(pancreas.combined, ndims = 50)

pancreas.combined <- FindNeighbors(pancreas.combined, reduction = "pca", dims = 1:20)
pancreas.combined <- FindClusters(pancreas.combined, resolution = 0.5)
pancreas.combined <- RunUMAP(pancreas.combined, reduction = "pca", 
                             dims = 1:20, min.dist = 0.5, n.neighbours = 100)

PCA_Seurat <- pancreas.combined[["pca"]]@cell.embeddings
UMAP_Seurat <- pancreas.combined[["umap"]]@cell.embeddings

#######
# ARI #
#######
ARI_Seurat <- adjustedRandIndex(metadata$CellType, pancreas.combined@meta.data$seurat_clusters)
write.table(ARI_Seurat, file = paste0("Output/ARI_",method,"_",proj,"_v",ver,".txt"), quote = FALSE)

#################
# Visualization #
#################
n.celltype <- length(unique(metadata$CellType))
celltype_names <- c("Alpha", "Beta", "Delta", "Gamma", "Acinar", "Ductal", "other")
# color_by_celltype <- c( "#fb8072","#fdb462","#C7EF85","#8dd3c7","#80b1d3","#bebada")
color_by_celltype_FACS <- brewer.pal(n.celltype, "Set3")

celltype_factor <- factor(metadata$CellType, levels = celltype_names)

### batch
B <- length(study_names)
color_by_batch <- viridis(B)
batch_factor <- factor(metadata$Study, levels = study_names)

### condition
condition_names <- c("ND", "T2D")
color_by_condition <- c("#ffffb3","#fbb4ae")
condition_factor <- factor(metadata$Disease, levels = condition_names)

df_lowD_Seurat <- data.frame(PC1 = PCA_Seurat[,1],
                      PC2 = PCA_Seurat[,2],
                      UMAP1 = UMAP_Seurat[,1], 
                      UMAP2 = UMAP_Seurat[,2], 
                      Batch = batch_factor,
                      CellType = celltype_factor,
                      Condition = condition_factor)



# Shuffle cells
set.seed(123)
cell_shuffle <- sample(1:N, N)
df_lowD_Seurat_shuffled <- df_lowD_Seurat[cell_shuffle, ]

###################
# Draw UMAP plots #
###################
type <- "UMAP"
feature1 <- "UMAP1"
feature2 <- "UMAP2"

## Color by batch
color_by <- "Batch"
color_group <- color_by_batch

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,".pdf"),width = 12, height = 8)
p <- ggplot(df_lowD_Seurat, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_manual(values = color_group) +
  xlab(feature1) + ylab(feature2) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333", 
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue", 
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2))
print(p)
dev.off()

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_shuffled.pdf"),width = 8, height = 8)
p <- ggplot(df_lowD_Seurat_shuffled, aes_string(x = feature1, y = feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_manual(values = color_group) +
  xlab(feature1) + ylab(feature2) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333", 
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue", 
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.position = "none",
    # legend.text = element_text(size = 24),
    # legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2))
print(p)
dev.off()

## Color by cell type
color_by <- "CellType"
color_group <- color_by_celltype_FACS
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,".pdf"),width = 12, height = 8)
p <- ggplot(df_lowD_Seurat, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_manual(values = color_group) +
  xlab(feature1) + ylab(feature2) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333", 
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue", 
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2))
print(p)
dev.off()

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_shuffled.pdf"),width = 8, height = 8)
p <- ggplot(df_lowD_Seurat_shuffled, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_manual(values = color_group) +
  xlab(feature1) + ylab(feature2) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333", 
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue", 
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.position = "none",
    # legend.text = element_text(size = 24),
    # legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2))
print(p)
dev.off()

## Color by cell type
color_by <- "Condition"
color_group <- color_by_condition

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,".pdf"),width = 12, height = 8)
p <- ggplot(df_lowD_Seurat, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_manual(values = color_group) +
  xlab(feature1) + ylab(feature2) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333", 
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue", 
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2))
print(p)
dev.off()

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_shuffled.pdf"),width = 8, height = 8)
p <- ggplot(df_lowD_Seurat_shuffled, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_manual(values = color_group) +
  xlab(feature1) + ylab(feature2) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333", 
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue", 
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.position = "none",
    # legend.text = element_text(size = 24),
    # legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2))
print(p)
dev.off()

#############################################
# Store the dimension reduction information #
#############################################
# Prepare the input for MELD
# Here, we only take the last two batches
PCA_output <- PCA_Seurat[metadata$Study %in% study_names[3:4], ]
UMAP_output <- UMAP_Seurat[metadata$Study %in% study_names[3:4], ]
metadata_output <-metadata[metadata$Study %in% study_names[3:4], ]

# rownames(PCA_output) <- paste0("Sample",metadata$sample, "_Cell", rep(1:ns[1], S))
colnames(PCA_output) <- paste0("PC",1:50)

# rownames(UMAP_output) <- rownames(PCA_output) 
colnames(UMAP_output) <- c("UMAP1", "UMAP2")

write.csv(PCA_output, file = paste0("RawCountData/PCA_Seurat_",proj,"_v",ver,".csv"), quote = FALSE)
write.csv(UMAP_output, file = paste0("RawCountData/UMAP_Seurat_",proj,"_v",ver,".csv"), quote = FALSE)
write.csv(metadata_output, file = paste0("RawCountData/metadata_",proj,"_v",ver,".csv"), quote = FALSE)

save(pancreas.combined, df_lowD_Seurat, color_by_batch, color_by_celltype_FACS, color_by_condition, 
     proj, method, ver,
     file = paste0("Output/",proj,"_low_dimension_",method,"_v",ver,".RData"))


##################
# MNN correction #
##################
method <- "MNN"

# Refer to https://bioconductor.org/packages/release/bioc/vignettes/batchelor/inst/doc/correction.html
# Create the separate singlecellExperiments object for each batch
sce <- list()
for(b in 1:B){
  
  cell_index <- metadata$Study == study_names[b]
  sce[[b]] <- SingleCellExperiment(assays = list(counts = y_obs[,cell_index],
                                                 logcounts = log1p(y_obs[,cell_index])),
                                   colData = metadata[cell_index,])
}

# Combine data sets from 6 batches together without correction
combined <- correctExperiments(Batch1 = sce[[1]], Batch2 = sce[[2]], 
                               Batch3 = sce[[3]], Batch4 = sce[[4]],
                               PARAM=NoCorrectParam())

# batch effect correction by fastMNN()
fastCorrect <- fastMNN(combined, batch = combined$batch)
fastCorrect <- runPCA(fastCorrect, dimred = "corrected", name = "pca.corrected")

# Generate UMAP
fastCorrect <- runUMAP(fastCorrect, dimred = "pca.corrected", pca = 20,
                       name = "UMAP", n_neighbors = 100)

PCA_MNN <- reducedDim(fastCorrect, type = "pca.corrected")
UMAP_MNN <- reducedDim(fastCorrect, type = "UMAP")

#################
# Visualization #
#################
n.celltype <- length(unique(metadata$CellType))
celltype_names <- c("Alpha", "Beta", "Delta", "Gamma", "Acinar", "Ductal", "other")
# color_by_celltype <- c( "#fb8072","#fdb462","#C7EF85","#8dd3c7","#80b1d3","#bebada")
color_by_celltype_FACS <- brewer.pal(n.celltype, "Set3")

celltype_factor <- factor(metadata$CellType, levels = celltype_names)

### batch
B <- length(study_names)
color_by_batch <- viridis(B)
batch_factor <- factor(metadata$Study, levels = study_names)

### condition
condition_names <- c("ND", "T2D")
color_by_condition <- c("#ffffb3","#fbb4ae")
condition_factor <- factor(metadata$Disease, levels = condition_names)

df_lowD_MNN <- data.frame(PC1 = PCA_MNN[,1],
                          PC2 = PCA_MNN[,2],
                          UMAP1 = UMAP_MNN[,1], 
                          UMAP2 = UMAP_MNN[,2], 
                          Batch = batch_factor,
                          CellType = celltype_factor,
                          Condition = condition_factor)

# Shuffle cells
set.seed(123)
cell_shuffle <- sample(1:N, N)
df_lowD_MNN_shuffled <- df_lowD_MNN[cell_shuffle, ]

###################
# Draw UMAP plots #
###################
type <- "UMAP"
feature1 <- "UMAP1"
feature2 <- "UMAP2"

## Color by batch
color_by <- "Batch"
color_group <- color_by_batch

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,".pdf"),width = 12, height = 8)
p <- ggplot(df_lowD_MNN, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_manual(values = color_group) +
  xlab(feature1) + ylab(feature2) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333", 
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue", 
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2))
print(p)
dev.off()

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_shuffled.pdf"),width = 8, height = 8)
p <- ggplot(df_lowD_MNN_shuffled, aes_string(x = feature1, y = feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_manual(values = color_group) +
  xlab(feature1) + ylab(feature2) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333", 
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue", 
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.position = "none",
    # legend.text = element_text(size = 24),
    # legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2))
print(p)
dev.off()

## Color by cell type
color_by <- "CellType"
color_group <- color_by_celltype_FACS
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,".pdf"),width = 12, height = 8)
p <- ggplot(df_lowD_MNN, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_manual(values = color_group) +
  xlab(feature1) + ylab(feature2) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333", 
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue", 
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2))
print(p)
dev.off()

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_shuffled.pdf"),width = 8, height = 8)
p <- ggplot(df_lowD_MNN_shuffled, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_manual(values = color_group) +
  xlab(feature1) + ylab(feature2) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333", 
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue", 
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.position = "none",
    # legend.text = element_text(size = 24),
    # legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2))
print(p)
dev.off()

## Color by cell type
color_by <- "Condition"
color_group <- color_by_condition

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,".pdf"),width = 12, height = 8)
p <- ggplot(df_lowD_MNN, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_manual(values = color_group) +
  xlab(feature1) + ylab(feature2) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333", 
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue", 
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2))
print(p)
dev.off()

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_shuffled.pdf"),width = 8, height = 8)
p <- ggplot(df_lowD_MNN_shuffled, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_manual(values = color_group) +
  xlab(feature1) + ylab(feature2) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333", 
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue", 
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.position = "none",
    # legend.text = element_text(size = 24),
    # legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2))
print(p)
dev.off()

#############################################
# Store the dimension reduction information #
#############################################
PCA_output <- PCA_MNN[metadata$Study %in% study_names[3:4], ]
UMAP_output <- UMAP_MNN[metadata$Study %in% study_names[3:4], ]
metadata_output <-metadata[metadata$Study %in% study_names[3:4], ]

colnames(metadata_output) <- c("CellID", "Sample", "Batch", "Condition", "CellType")
metadata_output$Cluster <- as.numeric(celltype_factor)[metadata$Study %in% study_names[3:4]]


# rownames(PCA_output) <- paste0("Sample",metadata$sample, "_Cell", rep(1:ns[1], S))
colnames(PCA_output) <- paste0("PC",1:ncol(PCA_MNN))

# rownames(UMAP_output) <- rownames(PCA_output) 
colnames(UMAP_output) <- c("UMAP1", "UMAP2")

write.csv(PCA_output, file = paste0("RawCountData/PCA_",method,"_",proj,"_v",ver,".csv"), quote = FALSE)
write.csv(UMAP_output, file = paste0("RawCountData/UMAP_",method,"_",proj,"_v",ver,".csv"), quote = FALSE)
write.csv(metadata_output, file = paste0("RawCountData/metadata_",proj,"_v",ver,".csv"), quote = FALSE)

save(fastCorrect, df_lowD_MNN, color_by_batch, color_by_celltype_FACS, color_by_condition, 
     proj, method, ver,
     file = paste0("Output/",proj,"_low_dimension_",method,"_v",ver,".RData"))
