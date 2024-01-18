## Conduct DA analysis on the pancreas dataset
rm(list=ls())

library(viridis)
library(ggplot2)
proj <- "pancreas"
ver <- 1
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

#############
# Draw UMAP #
#############
color_by <- "Chd_likelihood"
type <- "UMAP"
feature1 <- "UMAP1"
feature2 <- "UMAP2"
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,".pdf"),width = 12, height = 8)
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
