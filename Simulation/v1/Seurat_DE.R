# Apply Seurat to the simulated datasets
library(Seurat)
library(mclust)
library(viridis)
library(ggplot2)
library(pROC)

proj <- "simulation"
ver <- 1

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

simulation.combined <- FindNeighbors(simulation.combined, reduction = "pca", dims = 1:10)
simulation.combined <- FindClusters(simulation.combined, resolution = 0.5)
simulation.combined <- RunUMAP(simulation.combined, dims = 1:10)

#######
# ARI #
#######
ARI <- adjustedRandIndex(simulation.combined@meta.data$celltype, simulation.combined@meta.data$seurat_clusters)
write.table(ARI, file = paste0("Output/ARI_",method,"_",proj,"_v",ver,".txt"), quote = FALSE)


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
  }
  print(paste0("Finish identify DE genes for cell type ", k + 1, "."))
}


DE_ROC_Seurat <- data.frame(DE = as.vector(ifelse(true_DE, "DE", "Not DE")),
                            pval = as.vector(padj_Seurat))

roc_DE_Seurat <- roc(DE ~ pval, DE_ROC_Seurat, smooth = TRUE,
                     plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=FALSE,
                     auc.polygon.col="skyblue", print.auc=TRUE, xlim = c(1,0))

write.table(DE_ROC_Seurat, file = paste0("Output/",method,"_cell_type_specific_DE_genes_for_ROC_v",ver,".txt"))
