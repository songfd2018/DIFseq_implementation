##################################
# Work on the simulation dataset #
##################################
library(batchelor)
library(SingleCellExperiment)
library(dplyr)
library(patchwork)
library(scater)
library(miloR)
# library(Seurat)
set.seed(1234)
proj <- "covid_sc"
ver <- 0

# dir <- "D:/CUHKSZ/BUTseq/code/simulation/v60/"
# dir <- "/home/ondemand/songfangda/simulation/"
dir <- paste0("/home/ondemand/songfangda/DIFseq/",proj,"_comparison/Milo")
setwd(dir)

########################################
# load the raw count data and metadata #
########################################
# Load raw count data
y_obs <- read.table(paste0("../RawCountData/count_data_",proj,"_v",ver,".txt"))
y_obs <- as.matrix(y_obs)

# Load dimension
dim <- read.table(paste0("../RawCountData/dim_",proj,"_v",ver,".txt"))
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
metadata <- read.table(paste0("../RawCountData/metadata_",proj,"_v",ver,".txt"), header = TRUE)
# colnames(metadata) <- c("Cell","Donor","Batch", "Treatment", "CellType")

# metadata$CellType <- factor(metadata$CellType, levels = c("Alpha", "Beta", "Delta", "Gamma", "Acinar", "Ductal", "other"))
# metadata$Study <- factor(metadata$Study, levels = unique(metadata$Study))
# metadata$Disease <- factor(metadata$Disease, levels = unique(metadata$Disease))

# Load the gene list
gene_list <- unlist(read.table(paste0("../RawCountData/gene_list_",proj,"_v",ver,".txt"),stringsAsFactors = F))

# Add row and column names for the raw count data matrix
rownames(y_obs) <- gene_list
colnames(y_obs) <- metadata$SampleID

# indices of cells
btp_ind <- read.table(paste0("../RawCountData/tbinfor_",proj,"_v",ver,".txt"))

metadata$Pair <- btp_ind[,3]

b_ind <- btp_ind[,1]
nb <- table(b_ind)
cell_index <- NULL
for(b in 1:B){
  cell_index <- c(cell_index, 1:nb[b])
}

##################################
# batch effect correction by MNN #
##################################
# Refer to https://bioconductor.org/packages/release/bioc/vignettes/batchelor/inst/doc/correction.html
# Create the separate singlecellExperiments object for each batch
sce <- list()
for(b in 1:B){
  sce[[b]] <- SingleCellExperiment(assays = list(counts = y_obs[,metadata$Batch==b],
                                                 logcounts = log1p(y_obs[,metadata$Batch==b])),
                                   colData = metadata[metadata$Batch==b,])
}

# Combine data sets from 6 batches together without correction
combined <- correctExperiments(Batch1 = sce[[1]], Batch2 = sce[[2]],
                               PARAM=NoCorrectParam())
# set.seed(1234)
# combined <- runPCA(combined)
# combined <- runTSNE(combined, dimred="PCA")
# plotTSNE(combined, colour_by="batch")


# batch effect correction by fastMNN()
fastCorrect <- fastMNN(combined, batch = combined$batch)
fastCorrect <- runPCA(fastCorrect, dimred = "corrected", name = "pca.corrected")

# Generate UMAP
fastCorrect <- runUMAP(fastCorrect, dimred = "pca.corrected", pca = 15,
                       name = "UMAP", n_neighbors = 100)
plotUMAP(fastCorrect, colour_by="batch")
plotUMAP(fastCorrect, colour_by="celltype")


UMAP_MiloR <- reducedDim(fastCorrect, type = "UMAP")
save(UMAP_MiloR, file = paste0("04UMAP_coord_MiloR_",proj,".RData"))

##################################
# Differential abundance testing #
##################################
# Refer to https://rawcdn.githack.com/MarioniLab/miloR/7c7f906b94a73e62e36e095ddb3e3567b414144e/vignettes/milo_gastrulation.html#1_Load_data

# Create a Milo object
covid_milo <- Milo(fastCorrect)
covid_milo

# Construct KNN graph
covid_milo <- buildGraph(covid_milo, k = 30, d = 30, reduced.dim = "pca.corrected")
covid_milo <- makeNhoods(covid_milo, prop = 0.1, k = 30, d=30, refined = TRUE, reduced_dims = "pca.corrected")
plotNhoodSizeHist(covid_milo)

# Counting cells in neighbourhoods
covid_milo <- countCells(covid_milo, meta.data = metadata, samples="Sample_id")
head(nhoodCounts(covid_milo))

sim_design <- metadata[,c("pair", "treatment", "batch")]

# # Take the most abundance cluster identity among neighborhood cells
# neigh_mat <- nhoods(sim_milo)
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
sim_design$batch <- as.factor(sim_design$batch) 
sim_design$treatment <- paste0("Treatment",sim_design$treatment)
sim_design$treatment <- as.factor(sim_design$treatment) 
sim_design <- distinct(sim_design)
rownames(sim_design) <- sim_design$pair

# 3.6 computing neighbourhood connectivity
sim_milo <- calcNhoodDistance(sim_milo, d=30, reduced.dim = "pca.corrected")

# 3.7 testing
da_results <- testNhoods(sim_milo, design = ~ treatment, design.df = sim_design, reduced.dim = "UMAP")
head(da_results)

da_results %>%
  arrange(SpatialFDR) %>%
  head() 

# save.image("miloR_simulation_part2.RData")

# 5.1 Automatic grouping of neighbourhoods
sim_milo <- buildNhoodGraph(sim_milo)

plotUMAP(sim_milo) + plotNhoodGraphDA(sim_milo, da_results, alpha=0.05) +
  plot_layout(guides="collect")

## Find groups
da_results <- groupNhoods(sim_milo, da_results, da.fdr = 0.35, max.lfc.delta = 10)
head(da_results)
table(da_results$NhoodGroup)

plotNhoodGroups(sim_milo, da_results, layout="UMAP") 

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

adjustedRandIndex(cluster_milo,metadata$celltype)

save(fastCorrect, file = "04Output_MNN_simulation_v72.RData")
save(sim_milo, da_results, cluster_milo, file = "04Output_MiloR_simulation_v72.RData")

# # 4 inspecting DA testing results
# ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
# 
# ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
#   geom_point() +
#   geom_hline(yintercept = 1)
# 
# 
# sim_milo <- buildNhoodGraph(sim_milo)
# 
# ## Plot single-cell tSNE
# tsne_pl <- plotReducedDim(sim_milo, dimred = "TSNE", colour_by="treatment", text_by = "celltype", 
#                           text_size = 3, point_size=0.5) +
#   guides(fill="none")
# 
# ## Plot neighbourhood graph
# nh_graph_pl <- plotNhoodGraphDA(sim_milo, da_results, layout="TSNE",alpha=0.1) 
# 
# tsne_pl + nh_graph_pl +
#   plot_layout(guides="collect")
# 
# # DA in certain cell types
# da_results <- annotateNhoods(sim_milo, da_results, coldata_col = "celltype")
# head(da_results)
# 
# plotDAbeeswarm(da_results, group.by = "celltype", alpha = 0.2)
# 
# # Finding markers of DA populations
# sim_milo <- logNormCounts(sim_milo)
# da_results$NhoodGroup <- as.numeric(da_results$SpatialFDR < 0.1 & da_results$logFC < 0)
# da_nhood_markers <- findNhoodGroupMarkers(sim_milo, da_results, subset.row = rownames(sim_milo)[1:10])