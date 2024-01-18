# Apply Seurat and MNN to correct batch effects for the simulated datasets
library(Seurat)
library(mclust)
library(viridis)
library(ggplot2)
library(pROC)

proj <- "simulation"
ver <- 2

#####################
# Seurat correction #
#####################
method <- "Seurat"
setwd(paste0("Simulation/v",ver,""))

# Load data
load(paste0(proj, "_countdata_v",ver,".RData"))

simulation.list <- NULL
study_names <- unique(metadata$batch)
nb <- table(factor(metadata$batch))
gene_list <- paste0("Genes-",1:G)

colnames(y) <- paste0("Cells-",1:N)
rownames(y) <- gene_list
rownames(metadata) <- paste0("Cells-",1:N)

cell_index <- 0
for(b in 1:B){
  simulation.list[[b]] <- CreateSeuratObject(y[,cell_index + 1:nb[b]], 
                                             project = paste("Batch",b), 
                                             meta.data = metadata[cell_index + 1:nb[b],])
  cell_index <- cell_index + nb[b]
}

simulation.list <- lapply(X = simulation.list, FUN = function(x) {
  x <- NormalizeData(x)
})

# Here, we set the features of interest as the gene set analyzed in DIFseq
simulation.anchors <- FindIntegrationAnchors(object.list = simulation.list, anchor.features = gene_list)
simulation.combined <- IntegrateData(anchorset = simulation.anchors)
DefaultAssay(simulation.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
simulation.combined <- ScaleData(simulation.combined, verbose = FALSE)
simulation.combined <- RunPCA(simulation.combined, npcs = 50, verbose = FALSE)

ElbowPlot(simulation.combined, ndims = 50)

simulation.combined <- FindNeighbors(simulation.combined, reduction = "pca", dims = 1:15)
simulation.combined <- FindClusters(simulation.combined, resolution = 0.5)
simulation.combined <- RunUMAP(simulation.combined, dims = 1:15)

#######
# ARI #
#######
ARI <- adjustedRandIndex(simulation.combined@meta.data$celltype, simulation.combined@meta.data$seurat_clusters)
write.table(ARI, file = paste0("Output/ARI_",method,"_",proj,"_v",ver,".txt"), quote = FALSE)

#################
# Visualization #
#################
PCA_Seurat <- simulation.combined[["pca"]]@cell.embeddings
UMAP_Seurat <- simulation.combined[["umap"]]@cell.embeddings

df_lowD_Seurat <- data.frame(UMAP1 = UMAP_Seurat[,1], 
                             UMAP2 = UMAP_Seurat[,2], 
                             Batch = paste("Batch",simulation.combined$batch),
                             CellType = paste("CellType",simulation.combined$celltype),
                             Condition = paste("Condition",simulation.combined$condition))

color_by_celltype<-c("#F88A7E", "#FFD87D", "#ABD978", "#8097D3", "#9C7ACE")
celltype_factor <- factor(w)

### batch
color_by_batch <- viridis(B)
batch_factor <- factor(b_infor)

### condition
# color_by_condition <- c( "#4C00FF80","#00FF4D80","#FFFF0080")
color_by_condition<-c("#fbb4ae","#b3cde3")
condition_factor <- factor(t_infor)

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
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_shuffled.pdf"),
    width = 8, height = 8)
p <- ggplot(df_lowD_Seurat_shuffled, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size = 0.5) + theme_classic() +
  scale_color_manual(values = color_group) +
  xlab(NULL) + ylab(NULL) +
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
color_group <- color_by_celltype
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_shuffled.pdf"),
    width = 8, height = 8)
p <- ggplot(df_lowD_Seurat_shuffled, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size = 0.5) + theme_classic() +
  scale_color_manual(values = color_group) +
  xlab(NULL) + ylab(NULL) +
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
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_shuffled.pdf"),
    width = 8, height = 8)
p <- ggplot(df_lowD_Seurat_shuffled, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size = 0.5) + theme_classic() +
  scale_color_manual(values = color_group) +
  xlab(NULL) + ylab(NULL) +
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

###############################
# Cell-type-specific DE genes #
###############################
# ROC curve for cell-type-specific DE genes
# E_syn <- matrix(0,G,Num_Cond * K)
E_syn <- eta.syn != 0
E_syn <- E_syn[, -(1:K)]

table(paste0("Type_",metadata$celltype,"_Cond",metadata$condition), simulation.combined@meta.data$seurat_clusters)

# Define cluster-condition pairs to identify cell-type-specific DE genes
label_Seurat <- simulation.combined@meta.data$seurat_clusters
simulation.combined$celltype.condition <- paste0(label_Seurat, "_Cond" ,simulation.combined$condition)
Idents(simulation.combined) <- "celltype.condition"
cc_names <- names(table((simulation.combined$celltype.condition)))

label_Seurat <- simulation.combined@meta.data$seurat_clusters
celltype_matching <- unlist(mapClass(label_Seurat,metadata$celltype)$aTOb)
Seurat_ntype <- length(unique(label_Seurat))

# For the first cell type
true_DE <- NULL
padj_Seurat <- NULL
type_condition_names <- NULL

# cell_names <- colnames(simulation.combined)
for (k in 0:Seurat_ntype) {
  t <- 2
  ref_cells <-
    simulation.combined$condition == 1 &
    simulation.combined$seurat_clusters == k
  com_cells <-
    simulation.combined$condition == t &
    simulation.combined$seurat_clusters == k
  if (sum(ref_cells) > 10 & sum(com_cells) > 10) {
    celltype.DE <- FindMarkers(
      simulation.combined,
      logfc.threshold = 0,
      ident.1 = paste0(k, "_Cond1"),
      ident.2 = paste0(k, "_Cond", t),
      verbose = FALSE
    )
    
    order_genes <-
      order(factor(rownames(celltype.DE), levels = paste0("Genes-", 1:G)))
    padj_Seurat <-
      cbind(padj_Seurat, celltype.DE$p_val_adj[order_genes])
    true_DE <-
      cbind(true_DE, E_syn[, celltype_matching[k + 1] + (t - 2) * K])
    
    type_condition_names <- c(
      type_condition_names,
      paste0("Seurat_", k, "_True_", celltype_matching[k +
                                                         1], "_Cond", t)
    )
    
  }
  print(paste0("Finish identify DE genes for cell type ", k + 1, "."))
}


DE_ROC_Seurat <- data.frame(DE = as.vector(ifelse(true_DE, "DE", "Not DE")),
                            pval = as.vector(padj_Seurat))

roc_DE_Seurat <- roc(DE ~ pval, DE_ROC_Seurat, smooth = TRUE,
                     plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=FALSE,
                     auc.polygon.col="skyblue", print.auc=TRUE, xlim = c(1,0))

write.table(DE_ROC_Seurat, file = paste0("Output/",method,"_cell_type_specific_DE_genes_for_ROC_v",ver,".txt"))

#############################################
# Store the dimension reduction information #
#############################################
PCA_output <- PCA_Seurat
UMAP_output <- UMAP_Seurat

rownames(PCA_output) <- paste0("Sample",metadata$sample, "_Cell", rep(1:ns[1], S))
colnames(PCA_output) <- paste0("PC",1:50)

rownames(UMAP_output) <- rownames(PCA_output) 
colnames(UMAP_output) <- c("UMAP1", "UMAP2")

write.csv(PCA_output, file = paste0("RawCountData/PCA_Seurat_",proj,"_v",ver,".csv"), quote = FALSE)
write.csv(UMAP_output, file = paste0("RawCountData/UMAP_Seurat_",proj,"_v",ver,".csv"), quote = FALSE)

save(simulation.combined, df_lowD_Seurat, color_by_batch, color_by_celltype, color_by_condition, 
     proj, method, ver,
     file = paste0("Output/",proj,"_low_dimension_",method,"_v",ver,".RData"))

##################
# MNN correction #
##################
rm(list=ls())
library(batchelor)
library(scater)

proj <- "simulation"
ver <- 2
method <- "MNN"
# Load data
load(paste0(proj, "_countdata_v",ver,".RData"))

gene_list <- paste0("Gene_",1:G)
colnames(y) <- paste0("Cells-",1:N)
rownames(y) <- gene_list
rownames(metadata) <- paste0("Cells-",1:N)

# Refer to https://bioconductor.org/packages/release/bioc/vignettes/batchelor/inst/doc/correction.html
# Create the separate singlecellExperiments object for each batch
sce <- list()
for(b in 1:B){
  sce[[b]] <- SingleCellExperiment(assays = list(counts = y[,metadata$batch==b],
                                                 logcounts = log1p(y[,metadata$batch==b])),
                                   colData = metadata[metadata$batch==b,])
}

# Combine data sets from 6 batches together without correction
combined <- correctExperiments(Batch1 = sce[[1]], Batch2 = sce[[2]], 
                               Batch3 = sce[[3]], Batch4 = sce[[4]],
                               PARAM=NoCorrectParam())

# batch effect correction by fastMNN()
fastCorrect <- fastMNN(combined, batch = combined$batch)
fastCorrect <- runPCA(fastCorrect, dimred = "corrected", name = "pca.corrected")

# Generate tSNE
# fastCorrect <- runTSNE(fastCorrect, dimred= "pca.corrected", pca = 15,
#                       name = "tSNE")

# Generate UMAP
fastCorrect <- runUMAP(fastCorrect, dimred = "pca.corrected", pca = 15,
                       name = "UMAP", n_neighbors = 10)

PCA_MNN <- reducedDim(fastCorrect, type = "pca.corrected")
UMAP_MNN <- reducedDim(fastCorrect, type = "UMAP")

df_lowD_MNN <- data.frame(PC1 = PCA_MNN[,1],
                          PC2 = PCA_MNN[,2],
                          UMAP1 = UMAP_MNN[,1], 
                          UMAP2 = UMAP_MNN[,2], 
                          Batch = paste("Batch",metadata$Batch),
                          CellType = paste("CellType",metadata$CellType),
                          Condition = paste("Condition",metadata$Condition))

#################
# Visualization #
#################
color_by_celltype<-c("#F88A7E", "#FFD87D", "#ABD978", "#8097D3", "#9C7ACE")
celltype_factor <- factor(w)

### batch
color_by_batch <- viridis(B)
batch_factor <- factor(b_infor)

### condition
color_by_condition<-c("#fbb4ae","#b3cde3")
condition_factor <- factor(t_infor)

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

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_shuffled.pdf"),
    width = 8, height = 8)
p <- ggplot(df_lowD_MNN_shuffled, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size = 0.5) + theme_classic() +
  scale_color_manual(values = color_group) +
  xlab(NULL) + ylab(NULL) +
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
color_group <- color_by_celltype

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_shuffled.pdf"),
    width = 8, height = 8)
p <- ggplot(df_lowD_MNN_shuffled, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size = 0.5) + theme_classic() +
  scale_color_manual(values = color_group) +
  xlab(NULL) + ylab(NULL) +
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

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_shuffled.pdf"),
    width = 8, height = 8)
p <- ggplot(df_lowD_MNN_shuffled, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size = 0.5) + theme_classic() +
  scale_color_manual(values = color_group) +
  xlab(NULL) + ylab(NULL) +
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
PCA_output <- PCA_MNN
UMAP_output <- UMAP_MNN

rownames(PCA_output) <- paste0("Sample",metadata$sample, "_Cell", rep(1:ns[1], S))
colnames(PCA_output) <- paste0("PC",1:ncol(PCA_MNN))

rownames(UMAP_output) <- rownames(PCA_output) 
colnames(UMAP_output) <- c("UMAP1", "UMAP2")

write.csv(PCA_output, file = paste0("RawCountData/PCA_",method,"_",proj,"_v",ver,".csv"), quote = FALSE)
write.csv(UMAP_output, file = paste0("RawCountData/UMAP_",method,"_",proj,"_v",ver,".csv"), quote = FALSE)

save(fastCorrect, df_lowD_MNN, color_by_batch, color_by_celltype, color_by_condition, 
     proj, method, ver,
     file = paste0("Output/",proj,"_low_dimension_",method,"_v",ver,".RData"))