rm(list=ls())
# heatmap.3 referring https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R
library(ggplot2)
proj <- "pancreas"
ver <- 161
method <- "MELD"
today <- format(Sys.Date(), "%m%d")
setwd("E:/Study/DIFseq/ModelComparison/0913Pancreas/MELD")

####################################
# 1. Correct raw count data by MNN #
####################################

#########################
# 2. Run MELD in python #
#########################

#################################
# 3. Analyze the output of MELD #
#################################
# Load metadata
metadata <- read.csv(paste0("metadata_MELD_",proj,"_v",ver,".csv"))

# Load UMAP coordinates
UMAP_output <- read.csv(paste0("../RawCountData/UMAP_MNN_",proj,"_v",ver,".csv"))

# Load PCs
PCA_output <- read.csv(paste0("../RawCountData/PCA_MNN_",proj,"_v",ver,".csv"))

metadata$PC1 <- PCA_output$PC1
metadata$PC2 <- PCA_output$PC2
metadata$UMAP1 <- UMAP_output$UMAP1
metadata$UMAP2 <- UMAP_output$UMAP2

# Load sample likelihood obtained from MELD
sample_like <- read.csv(paste0("sample_like_MELD_",proj,"_v",ver,".csv"))

# Color palette
library(viridis)
## cell type
color_by_celltype<-c("#F88A7E", "#FFD87D", "#ABD978", "#8097D3", "#9C7ACE")
celltype_factor <- factor(metadata$celltype)

## batch
B <- length(unique(metadata$batch))
color_by_batch <- viridis(B)
batch_factor <- factor(metadata$batch)

## treatment
color_by_condition<-c("#fbb4ae","#b3cde3")
treatment_factor <- factor(metadata$condition)

color_by <- "Chd_likelihood"

##################
# Draw PCA plots #
##################
type <- "PCA"
feature1 <- "PC1"
feature2 <- "PC2"
pdf(paste0("../Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,".pdf"),width = 12, height = 8)
p <- ggplot(metadata, aes_string(x = feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_gradient2(low = "#80FFFF", mid = "#e0e0e0",high = "#FF80FF", midpoint = 0.5) +
  xlab(feature1) + ylab(feature2) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333", 
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue", 
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2)) +
  labs(color = "Cond 1 Like")
print(p)
dev.off()

#############
# Draw UMAP #
#############
type <- "UMAP"
feature1 <- "UMAP1"
feature2 <- "UMAP2"
pdf(paste0("../Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,".pdf"),width = 12, height = 8)
p <- ggplot(metadata, aes_string(x = feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_gradient2(low = "#80FFFF", mid = "#e0e0e0",high = "#FF80FF", midpoint = 0.5) +
  xlab(feature1) + ylab(feature2) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333", 
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue", 
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2)) +
  labs(color = "Cond 1 Like")
print(p)
dev.off()

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

# Jitter plot
ggplot(metadata, aes(x= CellType, y= Chd_likelihood, color = Condition)) +
  geom_jitter(size = 0.5) + theme_classic() +
  ylim(0,1) +
  scale_color_manual(values = color_by_condition) +
  xlab("Cell Type") + ylab("ND Likelihood") +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333", 
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue", 
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2))


save.image(paste0(today,method,"_analysis_",proj,"_v",ver,".RData"))

