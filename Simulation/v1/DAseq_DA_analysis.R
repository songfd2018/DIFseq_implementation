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

rm(list=ls())
proj <- "simulation"
ver <- 1
method <- "DAseq"
setwd(paste0("Simulation/v",ver,""))

python2use <- "D:/miniconda3/"
GPU <- 5

#################################
# Run on the simulation dataset #
#################################
# Load metadata
metadata <- read.csv(paste0("RawCountData/metadata_",proj,"_v",ver,".csv"))

# Load UMAP coordinates
UMAP_embed <- read.csv(paste0("RawCountData/UMAP_corrected_",proj,"_v",ver,".csv"))
rownames(UMAP_embed) <- UMAP_embed$X
UMAP_embed <- UMAP_embed[,-1]
UMAP_embed <- as.matrix(UMAP_embed)

# Load PCs
PCA_embed <- read.csv(paste0("RawCountData/PCA_corrected_",proj,"_v",ver,".csv"))
rownames(PCA_embed) <- PCA_embed$X
PCA_embed <- PCA_embed[,-1]
PCA_embed <- as.matrix(PCA_embed)

# metadata$UMAP1 <- UMAP_embed[,1]
# metadata$UMAP2 <- UMAP_embed[,2]

##########################
# Differential abundance #
##########################
labels_cond1 <- metadata[metadata$Condition == "Condition_1", "Sample"]
labels_cond2 <- metadata[metadata$Condition == "Condition_2", "Sample"]

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

type <- "UMAP"
feature1 <- "UMAP1"
feature2 <- "UMAP2"
color_by <- "DA_measure"
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
  labs(color = "DA measure")
print(p)
dev.off()

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_with_legend.pdf"),
    width = 8, height = 8)
p <- ggplot(metadata, aes_string(x = feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_gradient2(low = "#f03b20", mid = "#e0e0e0",high = "#2c7fb8", midpoint = 0) +
  xlab(NULL) + ylab(NULL) +
  guides(color = guide_colourbar(barwidth = 20, barheight = 2)) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333", 
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue", 
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.position = "top",
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
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,".pdf"),width = 12, height = 6)
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
pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,".pdf"),width = 12, height = 6)
da_regions$da.region.plot
dev.off()

# Contingency table between true DA cells and identified DA cells
DA_cell_syn <- paste0(metadata$Condition, "_" , metadata$CellType)
table(DA_cell_syn, da_regions$da.region.label)
