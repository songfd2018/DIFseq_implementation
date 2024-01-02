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

# To check
# library(batchelor) # MNN correction
# library(patchwork)

rm(list=ls())
proj <- "pancreas"
ver <- 161
method <- "Milo"
today <- format(Sys.Date(), "%m%d")

setwd("E:/Study/DIFseq/ModelComparison/0913Pancreas/")

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

# plotReducedDim(sce_pancreas, colour_by="batch", dimred = "umap") 
# sce_pancreas <- runPCA(sce_pancreas)
# 
# head(reducedDim(sce_pancreas, "pca"))
# head(reducedDim(sce_pancreas, "PCA"))
# reducedDims(sce_pancreas)

# No need of batch effect correction for pancreas_v161
# ##################################
# # batch effect correction by MNN #
# ##################################
# # Refer to https://bioconductor.org/packages/release/bioc/vignettes/batchelor/inst/doc/correction.html
# # Create the separate singlecellExperiments object for each batch
# sce <- list()
# for(b in 1:B){
#   sce[[b]] <- SingleCellExperiment(assays = list(counts = y[,metadata$batch==b],
#                                                  logcounts = log(y[,metadata$batch==b] + 1)),
#                                    colData = metadata[metadata$batch==b,])
# }
# 
# # Combine data sets from 6 batches together without correction
# combined <- correctExperiments(Batch1 = sce[[1]], Batch2 = sce[[2]], Batch3 = sce[[3]], 
#                                Batch4 = sce[[4]], Batch5 = sce[[5]], Batch6 = sce[[6]],
#                                PARAM=NoCorrectParam())
# 
# # batch effect correction by fastMNN()
# fastCorrect <- fastMNN(combined, batch = combined$batch)
# fastCorrect <- runPCA(fastCorrect, dimred = "corrected", name = "pca.corrected")
# 
# # Generate UMAP
# fastCorrect <- runUMAP(fastCorrect, dimred = "pca.corrected", pca = 15,
#                        name = "UMAP", n_neighbors = 100)
# plotUMAP(fastCorrect, colour_by="batch")
# plotUMAP(fastCorrect, colour_by="celltype")
# 
# 
# UMAP_MiloR <- reducedDim(fastCorrect, type = "UMAP")
# save(UMAP_MiloR, file = paste0("Output/04UMAP_coord_MiloR_",proj,"_v",ver,".RData"))

##################################
# Differential abundance testing #
##################################
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
head(nhoodCounts(pancreas_milo))
dim(nhoodCounts(pancreas_milo))

# Defining experimental design
pancreas_design <- metadata_sub[,c("Pair", "Disease", "Study")]

# # Take the most abundance cluster identity among neighborhood cells
# neigh_mat <- nhoods(pancreas_milo)
# neigh_index <- summary(neigh_mat)$j
# neigh_cell <- summary(neigh_mat)$i
#
# # calculate an ARI
# N_neigh <- max(neigh_index)
# neigh_celltype_count <- matrix(0, N_neigh, 5)
# for(ind in 1:length(neigh_index)){
#   # true cell type of the cell 
#   j <- neigh_index[ind]
#   celltype <- metadata$celltype[neigh_cell[j]]
#   neigh_celltype_count[j, celltype] <- neigh_celltype_count[j, celltype] + 1
# }
# head(neigh_celltype_count)
# # unfair to calculate ARI


## Convert batch info from integer to factor
pancreas_design$Study <- factor(pancreas_design$Study, levels = study_names[3:4]) 
pancreas_design$Disease <- as.factor(pancreas_design$Disease) 
pancreas_design <- distinct(pancreas_design)
rownames(pancreas_design) <- pancreas_design$Pair

# 3.6 computing neighbourhood connectivity
pancreas_milo <- calcNhoodDistance(pancreas_milo, d=30, reduced.dim = "pca")

# 3.7 testing
da_results <- testNhoods(pancreas_milo, design = ~ Disease, design.df = pancreas_design, reduced.dim = "umap")
head(da_results)

da_results %>%
  arrange(SpatialFDR) %>%
  head() 

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


###############################################################
# UMAP for cells the in last two batches colored by cell type #
###############################################################
colnames(metadata_sub)[1:5] <- c("CellID", "Sample", "Batch", "Condition", "CellType")

metadata_sub$UMAP1 <- UMAP_embed[,1]
metadata_sub$UMAP2 <- UMAP_embed[,2]

n.celltype <- length(unique(metadata_sub$CellType))
celltype_names <- c("Alpha", "Beta", "Delta", "Gamma", "Acinar", "Ductal", "other")
color_by_celltype<- rainbow(n.celltype)
metadata_sub$CellType <- factor(metadata_sub$CellType, levels = celltype_names)
# celltype_factor <- factor(metadata$CellType, levels = celltype_names)

### batch
color_by_batch <- viridis(4)
metadata_sub$Batch <- factor(metadata_sub$Batch, levels = unique(metadata_sub$Batch))
# batch_factor <- factor(metadata$Study, levels = study_names)

### condition
# color_by_condition <- c( "#4C00FF80","#00FF4D80","#FFFF0080")
condition_names <- c("ND", "T2D")
color_by_condition<-c("#67a9cf","#ef8a62")
metadata_sub$Condition <- factor(metadata_sub$Condition, levels = condition_names)
# condition_factor <- factor(metadata$Disease, levels = condition_names)

type <- "UMAP_two_batches"
feature1 <- "UMAP1"
feature2 <- "UMAP2"

## Color by batch
color_by <- "Batch"
color_group <- color_by_batch
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,".pdf"),width = 12, height = 8)
p <- ggplot(metadata_sub, aes_string(x= feature1, y= feature2, color = color_by)) +
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
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,".pdf"),width = 12, height = 8)
p <- ggplot(metadata_sub, aes_string(x= feature1, y= feature2, color = color_by)) +
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

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,".pdf"),width = 12, height = 8)
p <- ggplot(metadata_sub, aes_string(x= feature1, y= feature2, color = color_by)) +
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

save.image(paste0("Output/",today,method,"_analysis_",proj,"_v",ver,".RData"))

####################################################
# Draw a beeswarm plot with known cell type labels #
####################################################
da_results <- annotateNhoods(pancreas_milo, da_results, coldata_col = "celltype")
head(da_results)

ggplot(da_results, aes(celltype_fraction)) + geom_histogram(bins=50)

da_results$celltype <- paste0("CellType",da_results$celltype)

type <- "beeswarm"
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,".pdf"),width = 12, height = 8)
plotDAbeeswarm(da_results, group.by = "celltype")
dev.off()

#######################################
# Automatic grouping of neighborhoods #
#######################################
## Find groups
da_results <- groupNhoods(pancreas_milo, da_results, da.fdr = 0.7, max.lfc.delta = 2)
head(da_results)

table(da_results$NhoodGroup)

type <- "UMAP_neighborhoods_grouping"
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,".pdf"),width = 12, height = 8)
plotNhoodGroups(pancreas_milo, da_results, layout="umap") 
dev.off()

type <- "beeswarm_grouping"
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,".pdf"),width = 12, height = 8)
plotDAbeeswarm(da_results, "NhoodGroup")
dev.off()

# Assign cells to each Nhood Group
matr_cell_nhood <- nhoods(pancreas_milo)
cell_index <- summary(matr_cell_nhood)$i
Nhood_index <- summary(matr_cell_nhood)$j
cluster_milo <- rep(NA, N)
for(i in 1:N){
  nhoods_i <- which(cell_index==i)
  if(length(nhoods_i)>0){
    group_temp <- da_results$NhoodGroup[Nhood_index[nhoods_i]]
    freq <- table(group_temp)
    cluster_milo[i] <- names(freq[which.max(freq)])
  }
}

adjustedRandIndex(cluster_milo,metadata$celltype)


########################################################
# Identify Differential expressed genes for each group #
########################################################
nhood_markers <- findNhoodGroupMarkers(pancreas_milo, da_results,  
                                       aggregate.samples = TRUE, sample_col = "sample")

head(nhood_markers)

# Do comparison with true DE genes
DE_genes <- list()
n_groups <- length(unique(da_results$NhoodGroup))

for(k in 1:n_groups){
  DE_index <- (abs(nhood_markers[,(k - 1) * 2 + 2]) > 1) & (nhood_markers[,(k - 1) * 2 + 3] < 0.01)
  DE_genes[[k]] <- nhood_markers$GeneID[DE_index]
}

DE_genes_milo <- unique(unlist(DE_genes))

gene_list <- paste0("Gene_",1:G)

D_index <- beta.syn != 0
true_DE_index <- apply(D_index, 1, sum) > 0
true_DE_genes <- gene_list[true_DE_index]

########################################
# Identify cell-type-specific DE genes #
########################################
table(metadata$celltype, cluster_milo)

# Group 3 and Group 6 in Milo corresponds to cell type 1
k <- 1
true_DE_eta_index <- gene_list[eta.syn[, k + K] != 0]

dge_3 <- testDiffExp(pancreas_milo, da_results, design = ~ condition, meta.data = data.frame(colData(pancreas_milo)),
                     # subset.row = rownames(pancreas_milo)[1:10], 
                     subset.nhoods=da_results$NhoodGroup=="3")

df_group3 <- dge_3$`3`
DE_type1_group3 <- rownames(df_group3)[df_group3$adj.P.Val < 0.01]
length(DE_type1_group3)

dge_6 <- testDiffExp(pancreas_milo, da_results, design = ~ condition, meta.data = data.frame(colData(pancreas_milo)),
                     # subset.row = rownames(pancreas_milo)[1:10],
                     subset.nhoods=da_results$NhoodGroup=="6")
df_group6 <- dge_6$`6`
DE_type1_group6 <- rownames(df_group6)[df_group6$adj.P.Val < 0.01]
length(DE_type1_group6)

type <- "Venn_DE_type1"
Venn_col <- brewer.pal(3, "Set3")
venn.diagram(
  x = list(true_DE_eta_index, DE_type1_group3, DE_type1_group6),
  category.names = c("Type 1 True" , "Milo Group 3", "Milo Group 6"),
  filename = paste0("Images/",type,"_",proj,"_",method,"_v",ver,".jpg"),
  
  lwd = 2,
  lty = 'blank',
  fill = Venn_col,
  
  output=TRUE
)

# Group 1 and Group 5 in Milo corresponds to cell type 1
k <- 2
true_DE_eta_index <-  gene_list[eta.syn[, k + K] != 0]

dge_1 <- testDiffExp(pancreas_milo, da_results, design = ~ condition, meta.data = data.frame(colData(pancreas_milo)),
                     # subset.row = rownames(pancreas_milo)[1:10],
                     subset.nhoods=da_results$NhoodGroup=="1")

df_group1 <- dge_1$`1`
DE_type2_group1 <- rownames(df_group1)[df_group1$adj.P.Val < 0.01]
length(DE_type2_group1)


dge_5 <- testDiffExp(pancreas_milo, da_results, design = ~ condition, meta.data = data.frame(colData(pancreas_milo)),
                     # subset.row = rownames(pancreas_milo)[1:10], 
                     subset.nhoods=da_results$NhoodGroup=="5")

df_group5 <- dge_5$`5`
DE_type2_group5 <- rownames(df_group5)[df_group3$adj.P.Val < 0.01]
length(DE_type2_group5)

type <- "Venn_DE_type2"
Venn_col <- brewer.pal(3, "Set3")
venn.diagram(
  x = list(true_DE_eta_index, DE_type2_group5),
  category.names = c("Type 2 True" , "Milo Group 5"),
  filename = paste0("Images/",type,"_",proj,"_",method,"_v",ver,".jpg"),
  
  lwd = 2,
  lty = 'blank',
  fill = Venn_col[1:2],
  
  output=TRUE
)

# Group 1 and Group 5 in Milo corresponds to cell type 1
k <- 3
true_DE_eta_index <-  gene_list[eta.syn[, k + K] != 0]
length(true_DE_eta_index)

dge_7 <- testDiffExp(pancreas_milo, da_results, design = ~ condition, meta.data = data.frame(colData(pancreas_milo)),
                     # subset.row = rownames(pancreas_milo)[1:10],
                     subset.nhoods=da_results$NhoodGroup=="7")

df_group7 <- dge_7$`7`
DE_type3_group7 <- rownames(df_group7)[df_group7$adj.P.Val < 0.01]
length(DE_type3_group7)


dge_8 <- testDiffExp(pancreas_milo, da_results, design = ~ condition, meta.data = data.frame(colData(pancreas_milo)),
                     # subset.row = rownames(pancreas_milo)[1:10], 
                     subset.nhoods=da_results$NhoodGroup=="8")

df_group8 <- dge_8$`8`
DE_type3_group8 <- rownames(df_group8)[df_group8$adj.P.Val < 0.01]
length(DE_type3_group8)

type <- "Venn_DE_type3"
Venn_col <- brewer.pal(3, "Set3")
venn.diagram(
  x = list(true_DE_eta_index, DE_type3_group7, DE_type3_group8),
  category.names = c("Type 3 True" , "Milo Group 7", "Milo Group 8"),
  filename = paste0("Images/",type,"_",proj,"_",method,"_v",ver,".jpg"),
  
  lwd = 2,
  lty = 'blank',
  fill = Venn_col,
  
  output=TRUE
)

# ROC curve for cell-type-specific DE genes
E_syn <- matrix(0,G,Num_Cond * K)
E_syn <- eta.syn != 0
E_syn <- E_syn[,1:((Num_Cond - 1) * K) + K]

# Combine the p-values of the same cell type
# Adjust the order of rows
nhood_markers$GeneID <- factor(nhood_markers$GeneID,levels = paste0("Gene_",1:G))
nhood_markers <- nhood_markers[order(nhood_markers$GeneID),] 

celltype_matching <- unlist(mapClass(cluster_milo,metadata$celltype)$aTOb)
K_est <- (ncol(nhood_markers) - 1)/2
p_adjust <- matrix(NA, G, K_est)
DE_true <- matrix(NA, G, K_est)
for(k in 1:K_est){
  p_adjust[,k] <- nhood_markers[, k * 2 + 1]
  DE_true[,k] <- E_syn[,celltype_matching[k]]
}

# Take the minimum value of the adjusted p values across neighborhood clusters for each cell type
DE_ROC_milo <- data.frame(DE = as.vector(ifelse(DE_true[,10], "DE", "Not DE")),
                          p.adj = as.vector(p_adjust[,10]))

method <- "Milo"
pdf(paste0("../Images/ROC_DE_",proj,"_",method,"_v",ver,".pdf"),width = 6, height = 6)
roc_DE_DIFseq <- roc(DE ~ p.adj, DE_ROC_milo, smooth = TRUE,
                     plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=FALSE,
                     auc.polygon.col="skyblue", print.auc=TRUE, xlim = c(1,0))
dev.off()

save.image(paste0("Output/",today,method,"_analysis_",proj,"_v",ver,"_FinishDE.RData"))


# # 4 inspecting DA testing results
# ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
# 
# ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
#   geom_point() +
#   geom_hline(yintercept = 1)
# 
# 
# pancreas_milo <- buildNhoodGraph(pancreas_milo)
# 
# ## Plot single-cell tSNE
# tsne_pl <- plotReducedDim(pancreas_milo, dimred = "TSNE", colour_by="condition", text_by = "celltype", 
#                           text_size = 3, point_size=0.5) +
#   guides(fill="none")
# 
# ## Plot neighbourhood graph
# nh_graph_pl <- plotNhoodGraphDA(pancreas_milo, da_results, layout="TSNE",alpha=0.1) 
# 
# tsne_pl + nh_graph_pl +
#   plot_layout(guides="collect")
# 
# # DA in certain cell types
# da_results <- annotateNhoods(pancreas_milo, da_results, coldata_col = "celltype")
# head(da_results)
# 
# plotDAbeeswarm(da_results, group.by = "celltype", alpha = 0.2)
# 
# # Finding markers of DA populations
# pancreas_milo <- logNormCounts(pancreas_milo)
# da_results$NhoodGroup <- as.numeric(da_results$SpatialFDR < 0.1 & da_results$logFC < 0)
# da_nhood_markers <- findNhoodGroupMarkers(pancreas_milo, da_results, subset.row = rownames(pancreas_milo)[1:10])