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
library(pROC)

rm(list=ls())
proj <- "simulation"
ver <- 3
method <- "Milo"

setwd(paste0("Simulation/v",ver,""))

#############
# Load data #
#############
load(paste0("simulation_countdata_v",ver,".RData"))

# Build a SingleCellExperiment object
sce_simulation <- SingleCellExperiment(assays = list(counts = y, logcounts = log1p(y)),
                                       colData = metadata)

# Load PCA and UMAP for visualization
PCA_embed <- read.csv(paste0("RawCountData/PCA_MNN_",proj,"_v",ver,".csv"))
UMAP_embed <- read.csv(paste0("RawCountData/UMAP_MNN_",proj,"_v",ver,".csv"))

rownames(PCA_embed) <- colnames(y)
rownames(UMAP_embed) <- colnames(y)

PCA_embed <- PCA_embed[,-1]
UMAP_embed <- UMAP_embed[,-1]

reducedDim(sce_simulation, "pca") <- as.matrix(PCA_embed)
reducedDim(sce_simulation, "umap") <- as.matrix(UMAP_embed)

##################################
# Differential abundance testing #
##################################
# Refer to https://rawcdn.githack.com/MarioniLab/miloR/7c7f906b94a73e62e36e095ddb3e3567b414144e/vignettes/milo_gastrulation.html#1_Load_data
# Create a Milo object
sim_milo <- Milo(sce_simulation)
sim_milo

# Construct KNN graph
sim_milo <- buildGraph(sim_milo, k = 30, d = 30, reduced.dim = "pca")
sim_milo <- makeNhoods(sim_milo, prop = 0.1, k = 30, d=30, refined = TRUE, reduced_dims = "pca")
plotNhoodSizeHist(sim_milo)

# Counting cells in neighbourhoods
sim_milo <- countCells(sim_milo, meta.data = metadata,samples="pair")
head(nhoodCounts(sim_milo))
dim(nhoodCounts(sim_milo))

# Defining experimental design
sim_design <- metadata[,c("pair", "condition", "batch")]

## Convert batch info from integer to factor
sim_design$batch <- as.factor(sim_design$batch)
sim_design$condition <- paste0("Condition",sim_design$condition)
sim_design$condition <- as.factor(sim_design$condition)
sim_design <- distinct(sim_design)
rownames(sim_design) <- sim_design$pair

# 3.6 computing neighbourhood connectivity
sim_milo <- calcNhoodDistance(sim_milo, d=30, reduced.dim = "pca")

# 3.7 testing

da_results_allcom <- testNhoods(sim_milo, design = ~ condition,
                               model.contrasts = c("conditionCondition2",
                                                   "conditionCondition3",
                                                   "conditionCondition2-conditionCondition3"),
                               design.df = sim_design, reduced.dim = "umap")


# head(da_results)

head(da_results_allcom)

colnames(da_results_allcom)[1:3] <- paste0("logFC_",c("Cond1_Cond2", "Cond1_Cond3", "Cond2_Cond3"))

# take the threshold as sptial FDR < 0.1
# sum(da_results$SpatialFDR < 0.1)

sum(da_results_allcom$SpatialFDR < 0.1)


#######################################################
# Draw UMAP of each individual cell via average logFC #
#######################################################

# Compute average logFC for each cell
matr_cell_nhood <- nhoods(sim_milo)
cell_index <- summary(matr_cell_nhood)$i
Nhood_index <- summary(matr_cell_nhood)$j
average_logFC <- matrix(NA, N, 3)
for(i in 1:N){
  nhoods_i <- which(cell_index==i)
  if(length(nhoods_i)>0){
    average_logFC[i,1] <- mean(da_results_allcom$logFC_Cond1_Cond2[Nhood_index[nhoods_i]])
    average_logFC[i,2] <- mean(da_results_allcom$logFC_Cond1_Cond3[Nhood_index[nhoods_i]])
    average_logFC[i,3] <- mean(da_results_allcom$logFC_Cond2_Cond3[Nhood_index[nhoods_i]])
  }
}

metadata$AvelogFC_12 <- average_logFC[,1]
metadata$AvelogFC_13 <- average_logFC[,2]
metadata$AvelogFC_23 <- average_logFC[,3]
metadata$UMAP1 <- UMAP_embed[,1]
metadata$UMAP2 <- UMAP_embed[,2]

# Shuffle cells
set.seed(123)
cell_shuffle <- sample(1:N, N)
metadata_shuffled <- metadata[cell_shuffle, ]

color_by_condition2 <- c("#31a354","#2c7fb8","#f03b20")

type <- "UMAP"
feature1 <- "UMAP1"
feature2 <- "UMAP2"
color_group <- color_by_condition2

color_by <- "AvelogFC_12"
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,".pdf"),
    width = 10, height = 8)
p <- ggplot(metadata, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_gradient2(low = color_group[2], mid = "#e0e0e0",high = color_group[1], midpoint = 0) +
  xlab(feature1) + ylab(feature2) +
  guides(color = guide_colourbar(barwidth = 20, barheight = 1)) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333",
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue",
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.position = "top",
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2)) +
  labs(color = "Ave. LFC")
print(p)
dev.off()

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_shuffled.pdf"),
    width = 8, height = 8)
p <- ggplot(metadata_shuffled, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_gradient2(low = color_group[2], mid = "#e0e0e0",high = color_group[1], midpoint = 0) +
  xlab(feature1) + ylab(feature2) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333",
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue",
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.position = "none",
    # legend.text = element_text(size = 24),
    # legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2)) +
  labs(color = "Ave. LFC")
print(p)
dev.off()

color_by <- "AvelogFC_13"
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,".pdf"),
    width = 10, height = 8)
p <- ggplot(metadata, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_gradient2(low = color_group[3], mid = "#e0e0e0",high = color_group[1], midpoint = 0) +
  xlab(feature1) + ylab(feature2) +
  guides(color = guide_colourbar(barwidth = 20, barheight = 1)) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333",
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue",
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.position = "top",
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2)) +
  labs(color = "Ave. LFC")
print(p)
dev.off()

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_shuffled.pdf"),
    width = 8, height = 8)
p <- ggplot(metadata_shuffled, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_gradient2(low = color_group[3], mid = "#e0e0e0",high = color_group[1], midpoint = 0) +
  xlab(feature1) + ylab(feature2) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333",
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue",
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.position = "none",
    # legend.text = element_text(size = 24),
    # legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2)) +
  labs(color = "Ave. LFC")
print(p)
dev.off()

color_by <- "AvelogFC_23"
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,".pdf"),
    width = 10, height = 8)
p <- ggplot(metadata, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_gradient2(low = color_group[3], mid = "#e0e0e0",high = color_group[2], midpoint = 0) +
  xlab(feature1) + ylab(feature2) +
  guides(color = guide_colourbar(barwidth = 20, barheight = 1)) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333",
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue",
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.position = "top",
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2)) +
  labs(color = "Ave. LFC")
print(p)
dev.off()

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_shuffled.pdf"),
    width = 8, height = 8)
p <- ggplot(metadata_shuffled, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_gradient2(low = color_group[3], mid = "#e0e0e0",high = color_group[2], midpoint = 0) +
  xlab(feature1) + ylab(feature2) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333",
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue",
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.position = "none",
    # legend.text = element_text(size = 24),
    # legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2)) +
  labs(color = "Ave. LFC")
print(p)
dev.off()

####################################################
# Draw a beeswarm plot with known cell type labels #
####################################################
da_results_allcom <- annotateNhoods(sim_milo, da_results_allcom, coldata_col = "celltype")
head(da_results_allcom)

ggplot(da_results_allcom, aes(celltype_fraction)) + geom_histogram(bins=50)

da_results_allcom$celltype <- paste0("CellType",da_results_allcom$celltype)

type <- "beeswarm"
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,".pdf"),width = 12, height = 8)
plotDAbeeswarm(da_results_allcom, group.by = "celltype")
dev.off()

#######################################
# Automatic grouping of neighborhoods #
#######################################
## Find groups
da_results_allcom$logFC <- da_results_allcom$logFC_Cond1_Cond2
da_results_allcom <- groupNhoods(sim_milo, da_results_allcom, da.fdr = 0.5,  max.lfc.delta = 2)
head(da_results_allcom)

table(da_results_allcom$NhoodGroup)

type <- "UMAP_neighborhoods_grouping"
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,".pdf"),width = 12, height = 8)
plotNhoodGroups(sim_milo, da_results, layout="umap")
dev.off()

type <- "beeswarm_grouping"
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,".pdf"),width = 12, height = 8)
plotDAbeeswarm(da_results, "NhoodGroup")
dev.off()

# Assign cells to each Nhood Group
matr_cell_nhood <- nhoods(sim_milo)
cell_index <- summary(matr_cell_nhood)$i
Nhood_index <- summary(matr_cell_nhood)$j
cluster_milo <- rep(NA, N)
for(i in 1:N){
  nhoods_i <- which(cell_index==i)
  if(length(nhoods_i)>0){
    group_temp <- da_results_allcom$NhoodGroup[Nhood_index[nhoods_i]]
    freq <- table(group_temp)
    cluster_milo[i] <- names(freq[which.max(freq)])
  }
}

ARI <- adjustedRandIndex(cluster_milo,metadata$celltype)
write.table(ARI, file = paste0("Output/ARI_",method,"_",proj,"_v",ver,".txt"), quote = FALSE)

# ROC curve for cell-type-specific DE genes
E_syn <- eta.syn != 0
E_syn <- E_syn[, -(1:K)]

gene_list <- paste0("Genes-",1:G)
rownames(sim_milo) <- gene_list

celltype_matching <- unlist(mapClass(cluster_milo,metadata$celltype)$aTOb)
celltype_matching <- celltype_matching[order(as.numeric(names(celltype_matching)))]

metadata$condition <- factor(metadata$condition)

K_est <- length(celltype_matching)
p_adjust <- NULL
DE_true <- NULL
type_condition_names <- NULL

for(k in 1:K_est){

  # Judge whether the condition exists in the cluster of a given neighborhood
  nhs <- nhoods(sim_milo)
  nhood.x <- da_results_allcom$NhoodGroup == k
  if(sum(nhood.x) > 1){
    nhood.gr.cells <- rowSums(nhs[, nhood.x, drop=FALSE]) > 0
  }else{
    nhood.gr.cells <- nhs[, nhood.x] > 0
  }
  meta.temp <- metadata[nhood.gr.cells, ]
  cell_number_cond <- table(meta.temp$condition)
  subset_neighbors <- da_results_allcom$NhoodGroup==k

  if(cell_number_cond[1] > 10 & cell_number_cond[2] > 10){
    dge <- testDiffExp(sim_milo, da_results_allcom, design = ~ condition,
                       meta.data = metadata,
                       model.contrasts = paste0("condition",2),
                       subset.nhoods = subset_neighbors)

    order_genes <- order(factor(rownames(dge[[1]]),levels = gene_list))
    p_adjust <- cbind(p_adjust, dge[[1]]$adj.P.Val[order_genes])
    matched_cell_type <- celltype_matching[k]
    DE_true <- cbind(DE_true, E_syn[,matched_cell_type])

    type_condition_names <- c(type_condition_names,
                              paste0("Milo_",k,"_True_",celltype_matching[k+1],"_Cond",t))
  }

  if(cell_number_cond[1] > 10 & cell_number_cond[3] > 10){
    dge <- testDiffExp(sim_milo, da_results_allcom, design = ~ condition,
                       meta.data = metadata,
                       model.contrasts = paste0("condition",3),
                       subset.nhoods = subset_neighbors)

    order_genes <- order(factor(rownames(dge[[1]]),levels = gene_list))
    p_adjust <- cbind(p_adjust, dge[[1]]$adj.P.Val[order_genes])
    matched_cell_type <- celltype_matching[k]
    DE_true <- cbind(DE_true, E_syn[,matched_cell_type + K])

    type_condition_names <- c(type_condition_names,
                              paste0("Milo_",k,"_True_",celltype_matching[k+1],"_Cond",t))
  }
}

rownames(p_adjust) <- gene_list
colnames(p_adjust) <- type_condition_names

# Take the minimum value of the adjusted p values across neighborhood clusters for each cell type
DE_ROC_milo <- data.frame(DE = as.vector(ifelse(DE_true, "DE", "Not DE")),
                          p.adj = as.vector(p_adjust))

roc_DE_DIFseq <- roc(DE ~ p.adj, DE_ROC_milo, smooth = TRUE,
                     plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=FALSE,
                     auc.polygon.col="skyblue", print.auc=TRUE, xlim = c(1,0))

write.table(DE_ROC_milo, paste0("Output/",method,"_cell_type_specific_DE_genes_for_ROC_v",ver,".txt"))

