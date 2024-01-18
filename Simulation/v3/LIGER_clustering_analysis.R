#################################
# Cell type clustering by LIGER #
#################################
library(rliger)
library(mclust)
library(viridis)
library(ggplot2)
library(bigmemory)

rm(list=ls())
proj <- "simulation"
ver <- 3
method <- "LIGER"
setwd(paste0("Simulation/v",ver,""))

#############
# Load data #
#############
y <- read.big.matrix(filename = paste0("RawCountData/count_data_",proj,"_v",ver,".txt"),
                     sep = " ", type = "integer")


dim <- read.table(paste0("RawCountData/dim_",proj,"_v",ver,".txt"))
dim <- unlist(dim)
N <- dim[1]
S <- dim[2]
G <- dim[3]
B <- dim[4]
Num_Cond <- dim[5]
P <- dim[6]

BT_pair <- matrix(dim[6 + 1:(3*P)], byrow = TRUE, nrow = P) 
colnames(BT_pair) <- c("Batch", "Condition", "n_bt")

# rownames(y) <- paste0("Gene_",1:G)

# Load metadata
metadata <- read.table(paste0("RawCountData/metadata_",proj,"_v",ver,".txt"))

exp_matr <- list()
for(b in 1:B){
  exp_matr[[b]] <- y[, metadata$batch == b]
  rownames(exp_matr[[b]]) <- paste0("Gene_",1:G)
}
names(exp_matr) <- paste0("Batch_",1:B)

simulation_liger <- createLiger(raw.data = exp_matr)
simulation_liger <- normalize(simulation_liger)
simulation_liger@var.genes <- paste0("Gene_",1:G)
simulation_liger <- scaleNotCenter(simulation_liger)

##############################
# Joint matrix factorization #
##############################
simulation_liger <- optimizeALS(simulation_liger, k = 20)

###############################################
# Quantile normalization and joint clustering #
###############################################
simulation_liger <- quantile_norm(simulation_liger, refine.knn = FALSE)
simulation_liger <- louvainCluster(simulation_liger, resolution = 0.25)

ARI <- adjustedRandIndex(simulation_liger@clusters, metadata$celltype)
write.table(ARI, file = paste0("Output/ARI_",method,"_",proj,"_v",ver,".txt"), quote = FALSE)