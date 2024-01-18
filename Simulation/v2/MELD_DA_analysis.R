rm(list=ls())
# heatmap.3 referring https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R
library(ggplot2)
proj <- "simulation"
ver <- 2
method <- "MELD"


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
metadata <- read.csv(paste0("Output/metadata_MELD_",proj,"_v",ver,".csv"))

# Load UMAP coordinates
UMAP_output <- read.csv(paste0("RawCountData/UMAP_MNN_",proj,"_v",ver,".csv"))

# Load PCs
PCA_output <- read.csv(paste0("RawCountData/PCA_MNN_",proj,"_v",ver,".csv"))

metadata$PC1 <- PCA_output$PC1
metadata$PC2 <- PCA_output$PC2
metadata$UMAP1 <- UMAP_output$UMAP1
metadata$UMAP2 <- UMAP_output$UMAP2

# Load sample likelihood obtained from MELD
sample_like <- read.csv(paste0("Output/sample_like_MELD_",proj,"_v",ver,".csv"))

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

#############
# Draw UMAP #
#############
method <- "MELD"
type <- "UMAP"
feature1 <- "UMAP1"
feature2 <- "UMAP2"
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,".pdf"),
    width = 8, height = 8)
p <- ggplot(metadata, aes_string(x = feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_gradient2(low = "#2c7fb8", mid = "#e0e0e0",high = "#f03b20", midpoint = 0.5) +
  xlab(NULL) + ylab(NULL) +
  guides(color = guide_colourbar(barwidth = 1, barheight = 20)) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333", 
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue", 
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.position = "none",
    # legend.text = element_text(size = 24),
    # legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2)) +
  labs(color = "Cond 1 Like")
print(p)
dev.off()

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_with_legend.pdf"),
    width = 8, height = 8)
p <- ggplot(metadata, aes_string(x = feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_gradient2(low = "#2c7fb8", mid = "#e0e0e0",high = "#f03b20", midpoint = 0.5) +
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
  labs(color = "Cond 1 Like")
print(p)
dev.off()


