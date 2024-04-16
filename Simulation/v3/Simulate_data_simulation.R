################################################
# Generate simulated data for model comparison #
################################################
rm(list=ls())
# heatmap.3 referring https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R
library(ggplot2)

# For drawing heatmap to visualize effect sizes
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(mclust)
library(cluster)
library(edgeR)
library(tidyr)
library(factoextra)

proj <- "simulation"
ver <- 3
set.seed(3241)
dir <- paste0("Simulation/v",ver)
setwd(dir)

if(!dir.exists("Images")){
  dir.create("Images")
}

if(!dir.exists("Output")){
  dir.create("Output")
}

######################################
# Set the Synthetic Parameter Values #
######################################
# The number of conditions
Num_Cond <- 3

# The number of batches under each condition
B <- 6

# The number of active pairs
P <- 10

# The number of donors in each pair
m <- 10 

# The number of donors
S <- P * m
cell_per_sample <- 1000

# The grid of observed cells
BT_pair <- data.frame(Batch = c(1,1,1,2,2,3,4,4,5,6),
                      Condition = c(1,2,3,1,2,1,2,3,2,3),
                      nbt = rep(cell_per_sample * m, P))

# The number of cells per batch
ns <- rep(cell_per_sample, S) 

b_infor <- rep(BT_pair[,1], BT_pair[,3])
t_infor <- rep(BT_pair[,2], BT_pair[,3])
p_infor <- rep(1:P, BT_pair[,3])
s_infor <- rep(1:S, ns)

# The total number of cells
N <- sum(ns)

# The number of genes
G <- 2000

# The number of cell types
K <- 5

#The first column of gamma.syn denotes the intercept of 
#the logistic regression for dropout events
#The second column of gamma.syn denotes the odds ratios 
#of the logistic regression for dropout events
gamma.syn<-matrix(0,B,2)
gamma.syn[1,]<-c(-0.5,-0.5)
gamma.syn[2,]<-c(-0.8,-0.5)
gamma.syn[3,]<-c(-0.5,-0.8)
gamma.syn[4,]<-c(-1,-0.5)
gamma.syn[5,]<-c(-1,-0.5)
gamma.syn[6,]<-c(-1,-0.2)

#the log-scale baseline expression levels
# effects.perturb <- round(rnorm(100,1,0.2),digits = 1)
alpha.syn<-rep(NA,G)
alpha.syn[1:(G/4)]<-rep(2,G/4)
alpha.syn[(G/4+1):(G/4*2)]<-rep(1,G/4)
alpha.syn[(G/4*2+1):(G/4*3)]<-rep(0.5,G/4)
alpha.syn[(G/4*3+1):G]<-rep(0,G/4)

alpha.syn[1:60] <- 4
alpha.syn[G/4 + 1:30] <- 3
alpha.syn[G/4 * 2 + 1:30] <- 2.5
alpha.syn[G/4 * 3 + 1:30] <- 2

#the cell-type effects 
beta.syn<-matrix(0,G,K)

#the first cell type is regarded as the reference cell type 
#without cell-type effects 
beta.syn[,1] <- 0

#the cell-type effects of the second cell type
beta.syn[1:60,2] <- -2
beta.syn[61:80, 2] <- 2
beta.syn[81:100, 2] <- 1

beta.syn[G/4 + 1:30,2] <- -2
beta.syn[G/4 + 31:50,2] <- 1
beta.syn[G/4 + 51:60,2] <- 2

beta.syn[G/4 * 2 + 1:30,2] <- -2
beta.syn[G/4 * 2 + 31:45,2] <- 1.5
beta.syn[G/4 * 2 + 46:60,2] <- 2.5

beta.syn[G/4 * 3 + 1:30,2] <- -2
beta.syn[G/4 * 3 + 31:40,2] <- 2
beta.syn[G/4 * 3 + 41:60,2] <- 2.2

#the cell-type effects of the third cell type
beta.syn[1:60,3] <- -2
beta.syn[101:110, 3] <- 2

beta.syn[G/4 + 1:30,3] <- -2
beta.syn[G/4 + 61:65,3] <- 2
beta.syn[G/4 + 66:90,3] <- 2.4

beta.syn[G/4 * 2 + 1:30,3] <- -2
beta.syn[G/4 * 2 + 61:75,3] <- 1.2
beta.syn[G/4 * 2 + 76:90,3] <- 2

beta.syn[G/4 * 3 + 1:30,3] <- -2
beta.syn[G/4 * 3 + 61:70,3] <- 1.6
beta.syn[G/4 * 3 + 71:90,3] <- 1.8

#the cell-type effects of the forth cell type
beta.syn[1:60, 4] <- -2
beta.syn[111:130, 4] <- 2

beta.syn[G/4 + 1:30,4] <- -2
beta.syn[G/4 + 91:105,4] <- 2.1
beta.syn[G/4 + 106:120,4] <- 2.3

beta.syn[G/4 * 2 + 1:30,4] <- -2
beta.syn[G/4 * 2 + 91:100,4] <- 2
beta.syn[G/4 * 2 + 101:120,4] <- 1.8

beta.syn[G/4 * 3 + 1:30,4] <- -2
beta.syn[G/4 * 3 + 91:110,4] <- 2
beta.syn[G/4 * 3 + 111:120,4] <- 1.4

#the cell-type effects of the fifth cell type
beta.syn[1:60,5] <- -2
beta.syn[131:150, 5] <- 2

beta.syn[G/4 + 1:30,5] <- -2
beta.syn[G/4 + 121:150,5] <- 1.2

beta.syn[G/4 * 2 + 1:30,5] <- -2
beta.syn[G/4 * 2 + 121:130,5] <- 2
beta.syn[G/4 * 2 + 131:150,5] <- 1

beta.syn[G/4 * 3 + 1:30,5] <- -2
beta.syn[G/4 * 3 + 121:125,5] <- 0.8
beta.syn[G/4 * 3 + 126:150,5] <- 1.8

# Index of intrinsic genes
D.syn <- apply(beta.syn != 0, 1, sum) > 0

# condition effects 
eta.syn<-matrix(0,G,Num_Cond * K)
magnitude_condition <- 1

# for the first cell type
# Both intrinsic and DE genes
eta.syn[1:5, K + 1] <- -magnitude_condition
# Only DE genes
eta.syn[151:160, K + 1] <- magnitude_condition * 1.2
# Both intrinsic and DE genes

eta.syn[G/4 + 1:10,K + 1] <- -magnitude_condition * 1.1
# Only DE genes
eta.syn[G/4 + 151:165, K + 1] <- magnitude_condition * 1.3

# Both intrinsic and DE genes
eta.syn[G/4 * 2 + 1:10,K + 1] <- -magnitude_condition 
sum(eta.syn[, K + 1] != 0)

# for the second cell type
# Both intrinsic and DE genes
eta.syn[61:65, K + 2] <- -magnitude_condition
# Only DE genes
eta.syn[151:165, K + 2] <- magnitude_condition

# Both intrinsic and DE genes
# Only DE genes
eta.syn[G/4 + 31:45,K + 2] <- magnitude_condition
eta.syn[G/4 + 166:170, K + 2] <- magnitude_condition * 1.1

# Only DE genes
# eta.syn[G/4 * 3 + 26:35,K + 2] <- magnitude_condition
sum(eta.syn[, K + 2] != 0)

# the third condition
# upregulated or downregulated highly expressed genes in the corresponding cell types
eta.syn[,K * 2 + 1:K] <- 0

# for the first cell type
# Both intrinsic and DE genes
eta.syn[6:10, 2 * K + 1] <- -magnitude_condition
eta.syn[G/4 + 11:15,2 * K + 1] <- magnitude_condition * 1.3
eta.syn[G/4 * 3 + 11:20,2 * K + 1] <- -magnitude_condition 

# Only DE genes
eta.syn[161:165, 2 * K + 1] <- -magnitude_condition * 0.8
eta.syn[G/4 + 166:170, 2 * K + 1] <- magnitude_condition * 1.1

sum(eta.syn[, 2 * K + 1] != 0)
table(eta.syn[, 2* K + 1] !=0, beta.syn[,1] != 0)

# for the second cell type
# Both intrinsic and DE genes
eta.syn[66:75, 2 * K + 2] <- -magnitude_condition * 1.3
eta.syn[G/4 + 46:55, 2 * K + 2] <- magnitude_condition * 1.5

# Only DE genes
eta.syn[171:180, 2 * K + 2] <- -magnitude_condition * 0.9
eta.syn[G/4 + 171:180, 2 * K + 2] <- magnitude_condition

sum(eta.syn[, 2 * K + 2] != 0)
table(eta.syn[, 2* K + 2] !=0, beta.syn[,2] != 0)

# for the third cell type
# Both intrinsic and DE genes
eta.syn[101:105, K * 2 + 3] <- -magnitude_condition
eta.syn[G/4 + 71:80,K * 2 + 3] <- magnitude_condition * 1.2
eta.syn[G/4 * 2 + 76:80,K * 2 + 3] <- -magnitude_condition
eta.syn[G/4 * 3 + 91:100,K * 2 + 3] <- magnitude_condition * 1.1

# Only DE genes
eta.syn[181:190,K * 2 + 3] <- -magnitude_condition * 1.4
eta.syn[G/4 + 181:190,K * 2 + 3] <- magnitude_condition * 0.9

sum(eta.syn[, K * 2 + 3] != 0)
table(eta.syn[, 2* K + 3] !=0, beta.syn[,3] != 0)

# Knowledge about control genes
control_genes <- rep(0, G)

#################
# batch effects #
nu.syn<-matrix(NA,G,B)

#the first batch is taken as the reference batch
#without batch effects
nu.syn[,1] <- 0

#the batch effect of the second batch
nu.syn[,2] <- rep(c(1,-1,2,-2),each = G/4)

#the batch effect of the third batch
nu.syn[,3] <- rep(c(-1,0,-2,-2),each = G/4)

#the batch effect of the forth batch
nu.syn[,4] <- rep(c(-1,-2,-3,1),each = G/4)

#the batch effect of the fifth batch
nu.syn[,5] <- rep(c(1,-1,-2,2),each = G/4)

#the batch effect of the sixth batch
nu.syn[,6] <- rep(c(3,2,1,1),each = G/4)


####################
# Cell size factor #

delta.syn <- rep(0, N)
delta.syn <- rep(NA, N)
nbt <- table(p_infor)

# the first cell in each batch-treatment pair is regarded as the reference cell
# with the cell-specific size factors being 0
delta.syn[1:(nbt[1]/2)] <- 0
delta.syn[(nbt[1]/2 + 1):nbt[1]] <- 0.3

# the second pair
delta.syn[nbt[1] + 1:(nbt[2]/4)] <- 0
delta.syn[nbt[1] + (nbt[2]/4+1):(nbt[2]/2)] <- 0.3
delta.syn[nbt[1] + (nbt[2]/2+1):nbt[2]] <- -0.2

# the third pair
delta.syn[sum(nbt[1:2]) + 1:(nbt[3]/4)] <- 0
delta.syn[sum(nbt[1:2]) + (nbt[3]/4+1):(nbt[3]/2)] <- -0.15
delta.syn[sum(nbt[1:2]) + (nbt[3]/2+1):nbt[3]] <- -0.30

# the forth pair
delta.syn[sum(nbt[1:3]) + 1:(nbt[4]/2)] <- 0
delta.syn[sum(nbt[1:3]) + (nbt[4]/2+1):nbt[4]] <- -0.2

# the fifth pair
delta.syn[sum(nbt[1:4]) + 1:(nbt[5]/4)] <- 0
delta.syn[sum(nbt[1:4]) + (nbt[5]/4+1):(nbt[5]/2)] <- -0.3
delta.syn[sum(nbt[1:4]) + (nbt[5]/2+1):nbt[5]] <- 0.2

# the sixth pair
delta.syn[sum(nbt[1:5]) + 1:(nbt[6]/4)] <- 0
delta.syn[sum(nbt[1:5]) + (nbt[6]/4+1):(nbt[6]/2)] <- 0.4
delta.syn[sum(nbt[1:5]) + (nbt[6]/2+1):nbt[6]] <- -0.05

# the seventh pair
delta.syn[sum(nbt[1:6]) + 1:(nbt[7]/2)] <- 0
delta.syn[sum(nbt[1:6]) + (nbt[7]/2+1):nbt[7]] <- 0.15

# the eighth pair
delta.syn[sum(nbt[1:7]) + 1:(nbt[8]/4)] <- 0
delta.syn[sum(nbt[1:7]) + (nbt[8]/4+1):(nbt[8]/2)] <- 0.3
delta.syn[sum(nbt[1:7]) + (nbt[8]/2+1):nbt[8]] <- 0.6

# the ninth pair
delta.syn[sum(nbt[1:8]) + 1:(nbt[9]/2)] <- 0
delta.syn[sum(nbt[1:8]) + (nbt[9]/2+1):nbt[9]] <- 0.4

# the ten pair
delta.syn[sum(nbt[1:9]) + 1:(nbt[10]/4)] <- 0
delta.syn[sum(nbt[1:9]) + (nbt[10]/4+1):(nbt[10]/2)] <- 0.4
delta.syn[sum(nbt[1:9]) + (nbt[10]/2+1):nbt[10]] <- -0.1

# batch-specific and gene-specific overdispersion parameters
phi.syn<-matrix(5, G, B) #mean 2 var 0.5
phi.syn[(G*0.4 + 1):(G * 0.8),] <- 3
phi.syn[(G*0.8 + 1):G,] <- 2

rdir <- function(n, xi){
  .K <- length(xi)
  res <- matrix(NA, n, .K)
  for(i in 1:n){
    gam <- rgamma(.K, xi)
    res[i,] <- gam / sum(gam)
  }
  return(res)
}

# the cell-type proportions in each batch
pi.syn <- matrix(NA, S, K)
xi.syn <- matrix(NA, P, K)
xi.syn[1, ] <- c(0.3, 0.2, 0.2, 0.2, 0.1)
xi.syn[2, ] <- c(0.3, 0.2, 0.3, 0.1, 0.1)
xi.syn[3, ] <- c(0.3, 0.05, 0.2, 0.35, 0.1)
xi.syn[4, ] <- c(0.3, 0.2, 0.2, 0.2, 0.1)
xi.syn[5, ] <- c(0.3, 0.2, 0.3, 0.1, 0.1)
xi.syn[6, ] <- c(0.3, 0.2, 0.2, 0.2, 0.1)
xi.syn[7, ] <- c(0.3, 0.2, 0.3, 0.1, 0.1)
xi.syn[8, ] <- c(0.3, 0.05, 0.2, 0.35, 0.1)
xi.syn[9, ] <- c(0.3, 0.2, 0.3, 0.1, 0.1)
xi.syn[10, ] <- c(0.3, 0.05, 0.2, 0.35, 0.1)


xi.syn <- xi.syn * 3000
s_ind <- 0
for(p in 1:P){
  pi.syn[s_ind + 1:m, ] <- rdir(m, xi.syn[p, ])
  s_ind <- s_ind + m
}

###############################################
# Simulate latent variables and Observed data #
###############################################
# the cell-type indicators of each cell
w <- NULL

for(s in 1:S){
  w <- c(w, apply(rmultinom(ns[s],1,pi.syn[s,]),2,function(x){which(x==1)}))
}

# the indicators for dropout events
z <- matrix(NA, G, N)

# the underlying true expression levels
x <- matrix(NA, G, N)

#the observed expression levels
y <- matrix(NA, G, N)

#the logarithm of mean expreesion level of each gene in each cell
log.mu <- matrix(NA, G, N)

#generate the latent variable and observed data
for(i in 1:N){
  
  # obtain the batch index
  b <- b_infor[i]
  t <- t_infor[i]
  p <- p_infor[i]
  s <- s_infor[i]
  
  log.mu[,i] <- alpha.syn + beta.syn[,w[i]] + eta.syn[,(t-1) * K + w[i]] + 
    nu.syn[,b] + delta.syn[i]
  
  x[,i]<-rnbinom(G, size = phi.syn[,b], mu = exp(log.mu[,i]))
  
  prob_drop <- exp(gamma.syn[b,1] + gamma.syn[b,2] * x[,i])/(1+exp(gamma.syn[b,1] + gamma.syn[b,2] * x[,i])) 
  z[,i] <- rbinom(G, size = 1, prob = prob_drop)
  
  y[,i] <- (1-z[,i]) * x[,i]
  
}

metadata <- data.frame(batch = b_infor,
                       condition = t_infor,
                       pair = p_infor,
                       sample = s_infor,
                       celltype = w)

if(!dir.exists("RawCountData")){
  dir.create("RawCountData")
}

options(scipen = 10)

# Prepare the input for C++ source code
write.table(y, file = paste0("RawCountData/count_data_",proj,"_v",ver,".txt"),row.names = FALSE,col.names = FALSE)

write.table(c(N,S,G,B,Num_Cond,P), file = paste0("RawCountData/dim_",proj,"_v",ver,".txt"),row.names = FALSE, col.names = FALSE,quote = FALSE)
write.table(BT_pair, file = paste0("RawCountData/dim_",proj,"_v",ver,".txt"),row.names = FALSE, col.names = FALSE, append = TRUE)

metadata <- data.frame(batch = b_infor,
                       condition = t_infor,
                       pair = p_infor,
                       sample = s_infor,
                       celltype = w)

write.table(metadata, file = paste0("RawCountData/metadata_",proj,"_v",ver,".txt"))

write.table(metadata[,1:4], file = paste0("RawCountData/tbinfor_",proj,"_v",ver,".txt"), row.names = FALSE, col.names = FALSE)

write.table(control_genes, file = paste0("RawCountData/control_genes_",proj,"_v",ver,".txt"), row.names = FALSE, col.names = FALSE)

#####################################################
# Initialize the cell type labels by robust k-means #
#####################################################
set.seed(1234)
K_sel <- 3:8
rownames(y) <- paste0("Gene_",1:G)
colnames(y) <- paste0("Cell_",1:N)

# normalization
y_norm <- cpm(y, log = TRUE)
obs_df <- data.frame(t(y_norm))

# scaling
obs_df <- scale(obs_df)

# clustering on the first batch
cells_in_batch1 <- which(metadata[,1] == 1)
random_cells <- sample(cells_in_batch1, 1000)
batch_df <- obs_df[random_cells, ]

for(K_initial in K_sel) {
  cluster_batch <- pam(batch_df, k = K_initial)
  
  # medoids
  Medoids_batch <- cluster_batch$medoids
  w_inital <- rep(NA, N)
  for (i in 1:N) {
    # L1_dist <- colSums(abs(t(Medoids_batch)- obs_df[i,]))
    Pearson_dist <- apply(Medoids_batch, 1, cor, obs_df[i, ])
    w_inital[i] <- which.max(Pearson_dist)
  }
  
  write.table(
    w_inital,
    file = paste0("RawCountData/w_initial_K", K_initial, "_",proj,"_v",ver,".txt"),
    row.names = FALSE,
    col.names = FALSE
  )
}

##########################################
# Prepare the input for MELD python code #
##########################################
rownames(y) <- paste0("Gene_",1:G)
colnames(y) <- paste0("Sample",metadata$sample, "_Cell", rep(1:ns[1], S))

metadata_df <- data.frame(Batch = paste0("Batch_",b_infor),
                       Condition = paste0("Condition_",t_infor),
                       Pair = paste0("Pair_",t_infor),
                       Sample = paste0("Sample_",s_infor),
                       CellType = paste0("Type_", w))
rownames(metadata_df) <- colnames(y)

write.csv(y, file = paste0("RawCountData/count_data_",proj,"_v",ver,".csv"), quote = FALSE)
write.csv(metadata_df, file = paste0("RawCountData/metadata_",proj,"_v",ver,".csv"), quote = FALSE)

################################
# Draw the UMAP for raw counts #
################################
method <- "Raw"

# Draw PCA and UMAP for corrected read counts
sim_Raw_obj <- CreateSeuratObject(counts = y, 
                                     meta.data = metadata, 
                                     project = "sim_Raw") 

# scaling
sim_Raw_obj <- ScaleData(sim_Raw_obj)

# Run PCA
sim_Raw_obj <- RunPCA(sim_Raw_obj, features = rownames(sim_Raw_obj))
ElbowPlot(sim_Raw_obj, ndims= 50)

# Run UMAP
sim_Raw_obj <- RunUMAP(sim_Raw_obj, dims = 1:15)

PCA_Raw <- sim_Raw_obj[["pca"]]@cell.embeddings
UMAP_Raw <- sim_Raw_obj[["umap"]]@cell.embeddings
metadata_df$UMAP1 <- UMAP_Raw[,1]
metadata_df$UMAP2 <- UMAP_Raw[,2]

# Shuffle cells
set.seed(123)
cell_shuffle <- sample(1:N, N)
metadata_shuffled <- metadata_df[cell_shuffle, ]

# Draw UMAP plots
type <- "UMAP"
feature1 <- "UMAP1"
feature2 <- "UMAP2"

## Color by batch
color_by <- "Batch"
color_by_batch <- viridis(B)
color_group <- color_by_batch

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_shuffled.pdf"),
    width = 8, height = 8)
p <- ggplot(metadata_shuffled, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_manual(values = color_group) +
  xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333", 
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue", 
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.position = "none",
    panel.background = element_rect(colour = "black",size = 2))
print(p)
dev.off()

## Color by cell type
color_by <- "CellType"
color_by_celltype <- c("#F88A7E", "#FFD87D", "#ABD978", "#8097D3", "#9C7ACE")
color_group <- color_by_celltype

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_shuffled.pdf"),
    width = 8, height = 8)
p <- ggplot(metadata_shuffled, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_manual(values = color_group) +
  xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333", 
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue", 
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.position = "none",
    panel.background = element_rect(colour = "black",size = 2))
print(p)
dev.off()

## Color by cell type
color_by <- "Condition"
color_by_condition <- c("#ccebc5","#b3cde3","#fbb4ae")
color_group <- color_by_condition

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_shuffled.pdf"),
    width = 8, height = 8)
p <- ggplot(metadata_shuffled, aes_string(x= feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_manual(values = color_group) +
  xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(# face = "bold", #color = "#993333", 
    size = 24), #angle = 45),
    axis.text.y = element_text(# face = "bold", #color = "blue", 
      size = 24),#, angle = 45))
    axis.title=element_text(size=32,face="bold"),
    legend.position = "none",
    panel.background = element_rect(colour = "black",size = 2))
print(p)
dev.off()

#########################################
# Store the dimension reduction results #
#########################################
PCA_output <- PCA_Raw[,1:30]
UMAP_output <- UMAP_Raw

rownames(PCA_output) <- paste0("Sample",metadata$sample, "_Cell", rep(1:ns[1], S))
colnames(PCA_output) <- paste0("PC",1:30)

rownames(UMAP_output) <- rownames(PCA_output) 
colnames(UMAP_output) <- c("UMAP1", "UMAP2")

rownames(metadata_df) <- rownames(PCA_output)

write.csv(PCA_output, file = paste0("RawCountData/PCA_Raw_",proj,"_v",ver,".csv"), quote = FALSE)
write.csv(UMAP_output, file = paste0("RawCountData/UMAP_Raw_",proj,"_v",ver,".csv"), quote = FALSE)
write.csv(metadata_df, file = paste0("RawCountData/metadata_",proj,"_v",ver,".csv"), quote = FALSE)


save.image(paste0(proj,"_countdata_v",ver,".RData"))