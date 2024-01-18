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
ver <- 2
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
da_results <- testNhoods(sim_milo, design = ~ condition, design.df = sim_design, reduced.dim = "umap")
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
sim_milo <- buildNhoodGraph(sim_milo)
nh_graph_pl <- plotNhoodGraphDA(sim_milo, da_results, layout="umap", alpha=0.1) 

type <- "UMAP_neighborhoods"
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,".pdf"),width = 12, height = 8)
nh_graph_pl
dev.off()

#######################################################
# Draw UMAP of each individual cell via average logFC #
#######################################################
# Compute average logFC for each cell
matr_cell_nhood <- nhoods(sim_milo)
cell_index <- summary(matr_cell_nhood)$i
Nhood_index <- summary(matr_cell_nhood)$j
average_logFC <- rep(NA, N)
for(i in 1:N){
  nhoods_i <- which(cell_index==i)
  if(length(nhoods_i)>0){
    
    average_logFC[i] <- mean(da_results$logFC[Nhood_index[nhoods_i]])
    # group_temp <- da_results$NhoodGroup[Nhood_index[nhoods_i]]
    # freq <- table(group_temp)
    # cluster_milo[i] <- names(freq[which.max(freq)])
  }
}

metadata$AvelogFC <- average_logFC
metadata$UMAP1 <- UMAP_embed[,1]
metadata$UMAP2 <- UMAP_embed[,2]

type <- "UMAP"
feature1 <- "UMAP1"
feature2 <- "UMAP2"
color_by <- "AvelogFC"
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,".pdf"),
    width = 8, height = 8)
p <- ggplot(metadata, aes_string(x = feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_gradient2(low = "#f03b20", mid = "#e0e0e0",high = "#2c7fb8", midpoint = 0) +
  xlab(NULL) + ylab(NULL) +
  guides(color = guide_colourbar(barwidth = 20, barheight = 1)) +
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

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_with_legend.pdf"),
    width = 8, height = 8)
p <- ggplot(metadata, aes_string(x = feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_gradient2(low = "#f03b20", mid = "#e0e0e0",high = "#2c7fb8", midpoint = 0) +
  xlab(NULL) + ylab(NULL) +
  guides(color = guide_colourbar(barwidth = 1, barheight = 20)) +
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

save.image(paste0("Output/",method,"_analysis_",proj,"_v",ver,".RData"))

####################################################
# Draw a beeswarm plot with known cell type labels #
####################################################
da_results <- annotateNhoods(sim_milo, da_results, coldata_col = "celltype")
head(da_results)

ggplot(da_results, aes(celltype_fraction)) + geom_histogram(bins=50)

da_results$celltype <- paste0("CellType",da_results$celltype)

#######################################
# Automatic grouping of neighborhoods #
#######################################
## Find groups
da_results <- groupNhoods(sim_milo, da_results, da.fdr = 0.7, max.lfc.delta = 2)
head(da_results)

table(da_results$NhoodGroup)

type <- "UMAP_neighborhoods_grouping"
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,".pdf"),width = 12, height = 8)
plotNhoodGroups(sim_milo, da_results, layout="umap") 
dev.off()


# Assign cells to each Nhood Group
matr_cell_nhood <- nhoods(sim_milo)
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

ARI <- adjustedRandIndex(cluster_milo,metadata$celltype)
write.table(ARI, file = paste0("Output/ARI_",method,"_",proj,"_v",ver,".txt"), quote = FALSE)

########################################
# Identify cell-type-specific DE genes #
########################################
table(metadata$celltype, cluster_milo)
# ROC curve for cell-type-specific DE genes
E_syn <- matrix(0,G,Num_Cond * K)
E_syn <- eta.syn != 0
E_syn <- E_syn[,-(1:K)]

gene_list <- paste0("Genes-",1:G)
rownames(sim_milo) <- gene_list

celltype_matching <- unlist(mapClass(cluster_milo,metadata$celltype)$aTOb)
celltype_matching <- celltype_matching[order(as.numeric(names(celltype_matching)))]

K_est <- length(celltype_matching)
p_adjust <- NULL
DE_true <- NULL

# rownames(p_adjust) <- gene_list
# colnames(p_adjust) <- celltype_matching
for(k in 1:K_est){
  
  subset_neighbors <- da_results$NhoodGroup==k
  
  if(sum(subset_neighbors) > 10){
    dge <- testDiffExp(sim_milo, da_results, design = ~ condition, meta.data = metadata,
                       # subset.row = rownames(sim_milo)[1:10],
                       subset.nhoods=subset_neighbors)
    
    order_genes <- order(factor(rownames(dge[[1]]),levels = gene_list))
    
    
    p_adjust <- cbind(p_adjust, dge[[1]]$adj.P.Val[order_genes])
    
    matched_cell_type <- celltype_matching[k]
    DE_true <- cbind(E_syn[,matched_cell_type], DE_true)
  }

}


# Take the minimum value of the adjusted p values across neighborhood clusters for each cell type
DE_ROC_milo <- data.frame(DE = as.vector(ifelse(DE_true, "DE", "Not DE")),
                          p.adj = as.vector(p_adjust))

roc_DE_milo <- roc(DE ~ p.adj, DE_ROC_milo, smooth = TRUE,
                     plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=FALSE,
                     auc.polygon.col="skyblue", print.auc=TRUE, xlim = c(1,0))

write.table(DE_ROC_milo, paste0("Output/",method,"_cell_type_specific_DE_genes_for_ROC_v",ver,".txt"))