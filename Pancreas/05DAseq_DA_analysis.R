############################################
# Differential abundance analysis by DAseq #
############################################
library(DAseq)
library(SingleCellExperiment)
library(scuttle)
library(ggplot2)
library(Seurat)
library(mclust)
library(pROC)
library(viridis)

rm(list=ls())
proj <- "pancreas"
ver <- 161
method <- "DAseq"

today <- format(Sys.Date(), "%m%d")

setwd("E:/Study/DIFseq/ModelComparison/0913Pancreas/DAseq/")

python2use <- "D:/miniconda3/"
GPU <- 5

#################################
# Run on the simulation dataset #
#################################
# Load metadata
metadata <- read.csv(paste0("../RawCountData/metadata_",proj,"_v",ver,".csv"))

# Load UMAP coordinates
UMAP_embed <- read.csv(paste0("../RawCountData/UMAP_Seurat_",proj,"_v",ver,".csv"))
rownames(UMAP_embed) <- UMAP_embed$X
UMAP_embed <- UMAP_embed[,-1]
UMAP_embed <- as.matrix(UMAP_embed)

# Load PCs
PCA_embed <- read.csv(paste0("../RawCountData/PCA_Seurat_",proj,"_v",ver,".csv"))
rownames(PCA_embed) <- PCA_embed$X
PCA_embed <- PCA_embed[,-1]
PCA_embed <- as.matrix(PCA_embed)

metadata$PC1 <- PCA_embed[,1]
metadata$PC2 <- PCA_embed[,2]
metadata$UMAP1 <- UMAP_embed[,1]
metadata$UMAP2 <- UMAP_embed[,2]

##########################
# Differential abundance #
##########################
labels_cond1 <- metadata[metadata$Condition == "ND", "Sample"]
labels_cond2 <- metadata[metadata$Condition == "T2D", "Sample"]

# Compute DA meature for all cells
da_cells <- getDAcells(
  X = PCA_embed,
  cell.labels = metadata$Sample,
  labels.1 = unique(labels_cond1),
  labels.2 = unique(labels_cond2),
  k.vector = seq(50, 500, 50),
  plot.embedding = UMAP_embed
)

# Draw the UMAP colored by DA measure
metadata$DA_measure <- da_cells$da.pred
range(da_cells$da.pred)

type <- "UMAP_DA_measure"
feature1 <- "UMAP1"
feature2 <- "UMAP2"
color_by <- "DA_measure"
  
pdf(paste0("../Images/",type,"_",proj,"_",method,"_v",ver,".pdf"),width = 12, height = 8)
p <- ggplot(metadata, aes_string(x= feature1, y= feature2, color = color_by)) +
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
  labs(color = "DA measure")
print(p)
dev.off()

# Identify DA cells by the threshold -0.8 and 0.8
da_cells <- updateDAcells(
  X = da_cells, pred.thres = c(-0.8,0.8),
  plot.embedding = UMAP_embed
)


# Draw the UMAP colored by DA subpopulations
type <- "UMAP_DA_subpop"
pdf(paste0("../Images/",type,"_",proj,"_",method,"_v",ver,".pdf"),width = 12, height = 8)
da_cells$da.cells.plot
dev.off()

# Get DA subpopulations
da_regions <- getDAregion(
  X = PCA_embed,
  da.cells = da_cells,
  cell.labels = metadata$Sample,
  labels.1 = unique(labels_cond1),
  labels.2 = unique(labels_cond2),
  resolution = 0.01,
  plot.embedding = UMAP_embed,
)

str(da_regions[1:2])
type <- "UMAP_DA_region"
pdf(paste0("../Images/",type,"_",proj,"_",method,"_v",ver,".pdf"),width = 12, height = 8)
da_regions$da.region.plot
dev.off()

# Contingency table between true DA cells and identified DA cells
DA_cell_syn <- paste0(metadata$Batch, "_" , metadata$CellType)
table(DA_cell_syn, da_regions$da.region.label)

table(metadata$Sample[da_regions$da.region.label > 0], da_regions$da.region.label[da_regions$da.region.label > 0])

table(paste0(metadata$Sample,"_",metadata$CellType)[da_regions$da.region.label > 0], da_regions$da.region.label[da_regions$da.region.label > 0])

###############################################################
# UMAP for cells the in last two batches colored by cell type #
###############################################################
n.celltype <- length(unique(metadata$CellType))
celltype_names <- c("Alpha", "Beta", "Delta", "Gamma", "Acinar", "Ductal", "other")
color_by_celltype<- rainbow(n.celltype)
metadata$CellType <- factor(metadata$CellType, levels = celltype_names)
# celltype_factor <- factor(metadata$CellType, levels = celltype_names)

### batch
color_by_batch <- viridis(4)
metadata$Batch <- factor(metadata$Batch, levels = unique(metadata$Batch))
# batch_factor <- factor(metadata$Study, levels = study_names)

### condition
# color_by_condition <- c( "#4C00FF80","#00FF4D80","#FFFF0080")
condition_names <- c("ND", "T2D")
color_by_condition<-c("#67a9cf","#ef8a62")
metadata$Condition <- factor(metadata$Condition, levels = condition_names)
# condition_factor <- factor(metadata$Disease, levels = condition_names)

type <- "UMAP_two_batches"
feature1 <- "UMAP1"
feature2 <- "UMAP2"

## Color by batch
color_by <- "Batch"
color_group <- color_by_batch
pdf(paste0("../Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,".pdf"),width = 12, height = 8)
p <- ggplot(metadata, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size = 0.5) + theme_classic() +
  scale_color_manual(values = color_group[3:4]) +
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

## Color by cell type
color_by <- "CellType"
color_group <- color_by_celltype
pdf(paste0("../Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,".pdf"),width = 12, height = 8)
p <- ggplot(metadata, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size = 0.5) + theme_classic() +
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

## Color by cell type
color_by <- "Condition"
color_group <- color_by_condition

pdf(paste0("../Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,".pdf"),width = 12, height = 8)
p <- ggplot(metadata, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size = 0.5) + theme_classic() +
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

save.image(paste0(today,method,"_analysis_",proj,"_v",ver,".RData"))

###########################
# Differential expression #
###########################
# Identify marker genes
y <- read.table(paste0("RawCountData/count_data_",proj,"_v",ver,".txt"))
y <- as.matrix(y)
colnames(y) <- metadata$X
rownames(y) <- paste0("Genes-",1:nrow(y))
logy <- log1p(y)

# This needs tensor flow
STG_markers <- STGmarkerFinder(
  X = logy,
  da.regions = da_regions,
  lambda = 1.5, n.runs = 5, return.model = T,
  python.use = python2use, GPU = GPU
)

# Draw Venn plot among true intrinsic genes, cell-type-specific DE genes and DE gene identified by DAseq
load(paste0(proj,"_countdata_v",ver,".RData"))

# True intrinsic genes
gene_list <- paste0("Genes-",1:G)
D_index <- beta.syn != 0
true_intrinsic_index <- apply(D_index, 1, sum) > 0
true_intrinsic_genes <- gene_list[true_intrinsic_index]

# True DE genes
E_index <- eta.syn[, K + 1:K] != 0
E_index <- as.data.frame(E_index)
true_DE_index <- lapply(E_index, which)
true_DE_genes <- list()
for(k in 1:3){
  true_DE_genes[[k]] <- gene_list[true_DE_index[[k]]]
}

table(DA_cell_syn, da_regions$da.region.label)

# Group 2 and Group 5 corresponds to cell type 1
DE_type_cur <- STG_markers$da.markers[["2"]]
DE_type_cur$adjusted_pval <- p.adjust(DE_type_cur$p_value, method = "BH")
DE_group2 <- rownames(subset(DE_type_cur, adjusted_pval > 0.01)) 

DE_type_cur <- STG_markers$da.markers[["5"]]
DE_type_cur$adjusted_pval <- p.adjust(DE_type_cur$p_value, method = "BH")
DE_group5 <- rownames(subset(DE_type_cur, adjusted_pval > 0.01)) 

type <- "Venn_celltype1"
Venn_col <- brewer.pal(4, "Set3")
venn.diagram(
  x = list(true_intrinsic_genes, true_DE_genes[[1]], DE_group2, DE_group5),
  category.names = c("True intrinsic" , "Ture DE type 1", "DE group 2", "DE group 5"),
  filename = paste0("Images/",type,"_",proj,"_",method,"_v",ver,".jpg"),
  
  lwd = 2,
  lty = 'blank',
  fill = Venn_col,
  
  output=TRUE
)

type <- "UMAP_pred_group5"
method <- "DAseq"

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,".pdf"),width = 12, height = 6)
plotCellScore(
  X = UMAP_embed, score = STG_markers$model[["5"]]$pred
)
dev.off()

# Group 1 and Group 6 corresponds to cell type 2
DE_type_cur <- STG_markers$da.markers[["1"]]
DE_type_cur$adjusted_pval <- p.adjust(DE_type_cur$p_value, method = "BH")
DE_group1 <- rownames(subset(DE_type_cur, adjusted_pval > 0.01)) 

DE_type_cur <- STG_markers$da.markers[["6"]]
DE_type_cur$adjusted_pval <- p.adjust(DE_type_cur$p_value, method = "BH")
DE_group6 <- rownames(subset(DE_type_cur, adjusted_pval > 0.01)) 

type <- "Venn_celltype2"
Venn_col <- brewer.pal(4, "Set3")
venn.diagram(
  x = list(true_intrinsic_genes, true_DE_genes[[2]], DE_group1, DE_group6),
  category.names = c("True intrinsic" , "Ture DE type 2", "DE group 1", "DE group 6"),
  filename = paste0("Images/",type,"_",proj,"_",method,"_v",ver,".jpg"),
  
  lwd = 2,
  lty = 'blank',
  fill = Venn_col,
  
  output=TRUE
)

# Group 3 and Group 4 corresponds to cell type 3
DE_type_cur <- STG_markers$da.markers[["3"]]
DE_type_cur$adjusted_pval <- p.adjust(DE_type_cur$p_value, method = "BH")
DE_group3 <- rownames(subset(DE_type_cur, adjusted_pval > 0.01)) 

DE_type_cur <- STG_markers$da.markers[["4"]]
DE_type_cur$adjusted_pval <- p.adjust(DE_type_cur$p_value, method = "BH")
DE_group4 <- rownames(subset(DE_type_cur, adjusted_pval > 0.01)) 

type <- "Venn_celltype3"
Venn_col <- brewer.pal(4, "Set3")
venn.diagram(
  x = list(true_intrinsic_genes, true_DE_genes[[3]], DE_group3, DE_group4),
  category.names = c("True intrinsic" , "Ture DE type 3", "DE group 3", "DE group 4"),
  filename = paste0("Images/",type,"_",proj,"_",method,"_v",ver,".jpg"),
  
  lwd = 2,
  lty = 'blank',
  fill = Venn_col,
  
  output=TRUE
)

table(DA_cell_syn, da_regions$da.region.label)

# ROC curve for cell-type-specific DE genes
E_syn <- matrix(0,G,Num_Cond * K)
E_syn <- eta.syn != 0
E_syn <- E_syn[,1:3 + K]
rownames(E_syn)


# Combine the p-values of the same cell type
alpha_values <- matrix(NA, G, 3)
celltype_matching <- unlist(mapClass(da_regions$da.region.label,metadata$celltype)$aTOb[-1])
list_matching <- list()
for(k in 1:3){
  list_matching[[k]] <- which(celltype_matching == k) 
}

head(STG_markers$da.markers$`1`$p_value)
head(STG_markers$model$`1`$alpha)

# Take the minimum value of the adjusted p values across neighborhood clusters for each cell type
for(k in 1:3){
  matched_nhood <- list_matching[[k]]
  temp_alpha <- NULL
  for(l in matched_nhood){
    alpha <- STG_markers$model[[l]]$alpha
    name_factor <- factor(names(alpha), levels = paste0("Genes-",1:G))
    alpha <- alpha[order(name_factor)]
    temp_alpha <- cbind(temp_alpha, alpha)
    print(head(temp_alpha))
  }
  alpha_values[,k] <- apply(temp_alpha, 1, max)
}

DE_ROC_DAseq <- data.frame(DE = as.vector(ifelse(E_syn, "DE", "Not DE")),
                           alpha = as.vector(alpha_values))

method <- "DAseq"
pdf(paste0("../Images/ROC_DE_",proj,"_",method,"_v",ver,".pdf"),width = 6, height = 6)
roc_DE_DIFseq <- roc(DE ~ alpha, DE_ROC_DAseq, 
                     plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=FALSE,
                     auc.polygon.col="skyblue", print.auc=TRUE, xlim = c(1,0))
dev.off()

save.image(paste0("Output/",today,method,"_analysis_",proj,"_v",ver,"_FinishDE.RData"))

