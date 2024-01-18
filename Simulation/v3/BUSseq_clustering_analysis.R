#################################
# Cell type clustering by LIGER #
#################################
library(BUSseq)
library(mclust)
library(viridis)
library(ggplot2)
library(bigmemory)

rm(list=ls())
proj <- "simulation"
ver <- 3

method <- "BUSseq"
setwd(paste0("Simulation/v",ver))

K <- 5
#############
# Load data #
#############
y_obs <- read.table(paste0("RawCountData/count_data_",proj,"_v",ver,".txt"))
metadata <- read.table(paste0("RawCountData/metadata_",proj,"_v",ver,".txt"))

B <- length(unique(metadata$batch))

count_data <- list()
for(b in 1:B){
  count_data[[b]] <- y_obs[, metadata$batch == b]
  
}
names(count_data) <- paste0("Batch_",1:B)

BUSseq_start <- Sys.time()
BUSseqfits_res <- BUSseq_MCMC(ObservedData = count_data, 

                              seed = 1234,
                              n.cores = 16,
                              n.celltypes = K,
                              n.iterations = 5000)

BUSseq_end <- Sys.time()

print("It takes ", difftime(BUSseq_end, BUSseq_start, units = "mins"), " to run BUSseq for simulation_v102.")

ARI <- adjustedRandIndex(celltypes(BUSseqfits_res), metadata$celltype)
write.table(ARI, file = paste0("Output/ARI_",method,"_",proj,"_v",ver,".txt"), quote = FALSE)

