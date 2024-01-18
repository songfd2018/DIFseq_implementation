#################################
# Cell type clustering by LIGER #
#################################
library(rliger)
library(mclust)
library(viridis)
library(ggplot2)

rm(list=ls())
proj <- "simulation"
ver <- 2

method <- "LIGER"
setwd(paste0("Simulation/v",ver,""))

#############
# Load data #
#############
load(paste0("simulation_countdata_v",ver,".RData"))

exp_matr <- list()
for(b in 1:B){
  exp_matr[[b]] <- y[, metadata$batch == b]
  
}
names(exp_matr) <- paste0("Batch_",1:B)

simulation_liger <- createLiger(raw.data = exp_matr)
simulation_liger <- normalize(simulation_liger)
simulation_liger@var.genes <- rownames(y)
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