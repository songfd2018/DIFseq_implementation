rm(list=ls())
# heatmap.3 referring https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R
library(mclust)
library(ggplot2)
library(edgeR)
library(tidyr)
library(viridis)
library(pheatmap)
library(RColorBrewer)
library(Seurat)
library(pROC)
library(xtable)

source("../heatmap3.R")

proj <- "simulation"
method <- "DIFseq"
ver <- 1
dir <- paste0("Simulation/v",ver)
setwd(dir)
load(paste0(proj,"_countdata_v",ver,".RData"))

##################################################
# Select the optimal number of cell types by BIC #
##################################################
k_sel <- 3:8
BIC.record <- rep(NA, length(k_sel))
for (k in k_sel) {
  # load BIC
  inference_dir <- paste0("Inference_K",k,"/")
  BIC <- unlist(read.table(paste0(inference_dir, "BIC.txt")))
  BIC.record[k - k_sel[1] + 1] <- BIC
}

# Draw the BIC plot
pdf(paste0("Images/BIC_plot_simulation_v",ver,".pdf"),width = 6, height = 8)
par(mar = c(5.1,6.1,4.1,2.1)) 
plot(k_sel,BIC.record,xlab= "K",ylab = "BIC",type="n",cex.axis=2,cex.lab=3)
points(k_sel,BIC.record,type="b",pch=19,cex=3)
dev.off()

# Consider the optimal number of cell types
K_opt <- k_sel[which.min(BIC.record)]
dir_inference <- paste0("Inference_K",K_opt,"/")

#############################
# Load parameter estimation #
#############################
# logistic parameters 
gamma_est <- read.table(paste0(dir_inference,"gamma_est.txt"))

# load cell type labels
w_est <- read.table(paste0(dir_inference,"w_est.txt"))
w_est <- unlist(w_est)
table(w_est, w)

# load alpha_est
alpha_est <- read.table(paste0(dir_inference,"alpha_est.txt"))
alpha_est <- unlist(alpha_est)

# load beta_est
beta_est <- read.table(paste0(dir_inference,"beta_est.txt"))
beta_est <- matrix(unlist(beta_est),G,K_opt)
logmu_est<-beta_est+alpha_est

# load eta_est
eta_est <- read.table(paste0(dir_inference,"eta_est.txt"))
eta_est <- matrix(unlist(eta_est),G,K_opt * Num_Cond)

# load nu_est #
nu_est <- read.table(paste0(dir_inference,"nu_est.txt"))
nu_est <- matrix(unlist(nu_est),G,B)

# load delta_est
delta_est <- read.table(paste0(dir_inference,"delta_est.txt"))
delta_est <- unlist(delta_est)

# load phi_est
phi_est <- read.table(paste0(dir_inference,"phi_est.txt"))
phi_est <- matrix(unlist(phi_est),G,B)

# load pi_est
pi_est <- read.table(paste0(dir_inference,"pi_est.txt"))
pi_est <- matrix(unlist(pi_est),S,K_opt)

# load p and tau0
ptau1 <- unlist(read.table(paste0(dir_inference,"ptau1_est.txt")))
tau0 <- 0.01
tau1_est <- ptau1[1]
pbeta_est <- ptau1[2]
peta_est <- ptau1[3]

###########################################################
# Switch the order of cell types to match the true labels #
###########################################################
auto_switch <- function(label_true, label_est, cluster_num){
  
  res <- rep(NA, cluster_num)
  names(res) <- 1:cluster_num
  
  table_label <- table(factor(label_true, levels = 1:cluster_num), factor(label_est, levels = 1:cluster_num))
  
  num_remain <- cluster_num
  
  while(num_remain > 1){
    max_cell <- which.max(table_label)
    max_row <- (max_cell - 1) %% num_remain + 1
    max_col <- (max_cell - 1) %/% num_remain + 1
    
    res[as.numeric(rownames(table_label)[max_row])] <- as.numeric(colnames(table_label)[max_col])
    
    if(num_remain == 2){
      res[as.numeric(rownames(table_label)[3 - max_row])] <- as.numeric(colnames(table_label)[3 - max_col])
    }else{
      table_label <- table_label[-max_row,]
      table_label <- table_label[,-max_col]
    }
    num_remain <- num_remain - 1
  }
  return(res)
}

ARI <- adjustedRandIndex(w, w_est)
write.table(ARI, file = paste0("Output/ARI_",method,"_",proj,"_v",ver,".txt"), quote = FALSE)

est_switch <- auto_switch(w, w_est, K_opt)
table_celltype <- table(w_est, metadata$celltype)
table_celltype <- table_celltype[est_switch,]
table_celltype

#####################################
# Barplot for cell type proportions #
#####################################
pi_ordered <- pi_est[,est_switch]

celltype_name <- paste0("Cell type ",1:K_opt)
sample_id <- paste0("Sample ",1:S)

colnames(pi_ordered) <- celltype_name
rownames(pi_ordered) <- sample_id

df_prop <- data.frame(
  Sample=factor(rep(sample_id, K_opt), levels = paste("Sample",1:100)),
  CellType = factor(rep(celltype_name, each = S)),
  Prop_est = as.vector(pi_ordered),
  Prop_true = as.vector(pi.syn))

color_by_celltype<-c("#F88A7E", "#FFD87D", "#ABD978", "#8097D3", "#9C7ACE")

pdf(file = paste0("Images/CellTypeProportion_syn_",proj,"_v",ver,".pdf"), width = 16, height = 5)
p <- ggplot(data = df_prop, mapping = aes(x = Sample, fill = CellType, y = Prop_true)) + 
  geom_col(width = 1, color = "#939393") +
  scale_fill_manual(values=color_by_celltype) +
  labs(x = "Sample", y = "Cell Type Proportion",colour = "Cell Type") +
  theme(axis.text=element_text(size=20), 
        axis.title=element_text(size=24,face="bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size=20, face = "bold"), #change legend title font size
        legend.text = element_text(size=20),
        panel.grid = element_blank(),
        panel.background =element_blank(),
        panel.border = element_blank())
p
dev.off()

pdf(file = paste0("Images/CellTypeProportion_est_",proj,"_v",ver,".pdf"), width = 16, height = 5)
p <- ggplot(data = df_prop, mapping = aes(x = Sample, fill = CellType, y = Prop_est)) + 
  geom_col(width = 1, color = "#939393") +
  scale_fill_manual(values=color_by_celltype) +
  labs(x = "Sample", y = "Cell Type Proportion",colour = "Cell Type") +
  theme(axis.text=element_text(size=20), 
        axis.title=element_text(size=24,face="bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size=20, face = "bold"), #change legend title font size
        legend.text = element_text(size=20),
        panel.grid = element_blank(),
        panel.background =element_blank(),
        panel.border = element_blank())
p
dev.off()

#################################
# Heatmap of biological effects #
#################################
scale01 <- function(x, low = min(x), high = max(x)) {
  x <- (x - low)/(high - low)
  x
}

draw_keys <- function(effects_matrix, color_bar, color_key = NULL, break_effects = NULL){
  if(is.null(color_key) & is.null(break_effects)){
    colorsChoice <- colorRampPalette(c("#F2F2F2","#060606"))
    color_key <- colorsChoice(10)
    len_breaks <- length(color_key) + 1
    break_effects <- seq(min(effects_matrix), max(effects_matrix), length.out = len_breaks)
  }else if(is.null(color_key)){
    colorsChoice <- colorRampPalette(c("#F2F2F2","#060606"))
    color_key <- colorsChoice(length(break_effects) - 1)
  }else if(is.null(break_effects)){
    len_breaks <- length(color_key) + 1
    break_effects <- seq(min(effects_matrix), max(effects_matrix), length.out = len_breaks)
  }

  z <- seq(min(break_effects), max(break_effects), length = length(color_key))
  image(z = matrix(z, ncol = 1), col = color_key, breaks = break_effects, xaxt = "n", yaxt = "n")
  lv <- pretty(break_effects)
  xv <- scale01(as.numeric(lv), min(break_effects), max(break_effects))
  axis(1, at = xv, labels = lv, cex.axis = 2.5)
  mtext(side = 1, "Value", line = 3.5, cex = 3)
  title(main = "Color Key", cex.main = 3)

}

heatmap_effects <- function(effects_matrix, color_bar, color_key = NULL, break_effects = NULL){

  if(is.null(color_key) & is.null(break_effects)){
    colorsChoice <- colorRampPalette(c("#F2F2F2","#060606"))
    color_key <- colorsChoice(10)
    len_breaks <- length(color_key) + 1
    break_effects <- seq(min(effects_matrix), max(effects_matrix), length.out = len_breaks)
  }else if(is.null(color_key)){
    colorsChoice <- colorRampPalette(c("#F2F2F2","#060606"))
    color_key <- colorsChoice(length(break_effects) - 1)
  }else if(is.null(break_effects)){
    len_breaks <- length(color_key) + 1
    break_effects <- seq(min(effects_matrix), max(effects_matrix), length.out = len_breaks)
  }

  out_fig <- heatmap.3(effects_matrix,
                       dendrogram = "none",#with cluster tree
                       Rowv = FALSE, Colv = FALSE,
                       labRow = FALSE, labCol = FALSE,
                       ColSideColors = color_bar,
                       lmat=rbind(c(5,4),c(0,1), c(3,2)),#1=heatmap, 2=row dendogram, 3=col dendogram, 4= key
                       lhei=c(0.3,0.4,3.6),
                       lwid = c(0.3,3),
                       col=color_key, breaks = break_effects, key = FALSE)

  return(out_fig)

}

# Heatmap for the log-scale mean expression levels---$\alpha_{g} + \beta_{gk}$
logmu_syn <- alpha.syn + beta.syn
logmu_est <- alpha_est + beta_est[, est_switch]

break_mu <- seq(-0.5,4.2, length.out = 101)

jpeg(paste0("Images/heatmap_logmu_syn_v",ver,".jpg"),width = 480, height = 960)
heat_logmu <- heatmap_effects(logmu_syn, 
                              color_by_celltype,
                              break_effects = break_mu)
dev.off()

jpeg(paste0("Images/colorkey_logmu_syn_v",ver,".jpg"),width = 480, height = 240)
draw_keys(mu.syn,
          color_by_celltype,
          break_effects = break_mu)
dev.off()

jpeg(paste0("Images/heatmap_logmu_est_v",ver,"_K",K_opt,".jpg"), width = 480, height = 960, quality = 100)
heat_logmu <- heatmap_effects(logmu_est, 
                              color_by_celltype, 
                              break_effects = break_mu)
dev.off()

# Heatmap for condition effects---$\eta_{tgk}$
color_by_condition<-c("#fbb4ae","#b3cde3")

eta_organized_syn <- eta.syn
eta_organized <- eta_est

for(k in 1:K_opt){
  for(t in 1:Num_Cond){
    
    eta_organized_syn[,(k-1) * Num_Cond + t] <- 
      eta.syn[,K_opt * (t-1) + k]
    
    eta_organized[,(k-1) * Num_Cond + t] <- 
      eta_est[,K_opt * (t-1) + est_switch[k]]
    
  }
}

# color bar
color_bar <- cbind(
  rep(color_by_condition,K_opt),
  rep(color_by_celltype,each = Num_Cond))

# color key
break_eta <- seq(-0.8,0.8, length.out = 101)
colorsChoice <- colorRampPalette(c("#1D78B2","#F2F2F2","#C4012D"))
treatment_key <- colorsChoice(100)


jpeg(paste0("Images/heatmap_eta_syn_v",ver,".jpg"),
     width = 480, height = 960)
heat_logmu <- heatmap_effects(eta_organized_syn, 
                              color_bar,
                              break_effects = break_eta, 
                              color_key = treatment_key)
dev.off()

jpeg(paste0("Images/colorkey_eta_syn_v",ver,".jpg"),width = 480, height = 240)
draw_keys(mu.syn,
          break_effects = break_eta, 
          color_key = treatment_key)
dev.off()

jpeg(paste0("Images/heatmap_eta_est_v",ver,".jpg"),
     width = 480, height = 960)
heat_logmu <- heatmap_effects(eta_organized, 
                              color_bar,
                              break_effects = break_eta, 
                              color_key = treatment_key)
dev.off()


######################################
# Generate the corrected read counts #
######################################
adjusted_values <- function(ReadCount, Indicators, .K,
                            .alpha, .beta, .eta, .nu, .delta, .phi, .w){
  CorrectedCount <- ReadCount
  N <- ncol(ReadCount)
  for(i in 1:N){
    
    b <- Indicators[i,1]
    t <- Indicators[i,2]
    w <- .w[i]
    
    # percentile
    px <- pnbinom(ReadCount[,i], size = .phi[,b], mu = exp(.alpha + .beta[,w] + .eta[,(t-1) * .K + w] + .nu[,b] + .delta[i]))
    pxminus1 <- pnbinom(ReadCount[,i] - 1, size = .phi[,b], mu = exp(.alpha + .beta[,w] + .eta[,(t-1) * .K + w] + .nu[,b] + .delta[i]))
    
    # get the aligned percentile
    local_u <- runif(G) * (px - pxminus1) + pxminus1
    local_u <- ifelse(local_u > 0.9999, 0.9999, local_u)
    
    # obtain the quantile
    CorrectedCount[,i] <- qnbinom(local_u, size = .phi[,1], mu = exp(.alpha + .beta[,w] + .eta[,(t-1) * .K + w]))
    
    if(i %% 1000 == 0){
      print(paste("Finish the correction of", i, "cells..."))
      
    }
  }
  return(CorrectedCount)
}

x_imputed <- read.table(paste0(dir_inference,"imputed_count.txt"))

start_time<-Sys.time()
message("Calculate corrected read counts:")

btp_ind <- cbind(b_infor,t_infor,p_infor)
x_corrected <-adjusted_values(x_imputed, btp_ind, K_opt,
                              alpha_est, beta_est, eta_est,nu_est,delta_est,phi_est,w_est)
write.table(x_corrected, file = paste0(dir_inference,"/corrected_count.txt"), row.names = FALSE, col.names = FALSE)
end_time<-Sys.time()
running_time<-difftime(end_time, start_time, units = "mins")
message("It takes ",running_time," mins to calculate the corrected data.")

colnames(x_corrected) <- paste0("Cells-",1:N)
rownames(x_corrected) <- paste0("Genes-",1:G)
rownames(metadata) <- paste0("Cells-",1:N)

########################################################
# Heatmap for raw count data and corrected read counts #
########################################################
set.seed(1234)
# Randomly sample 500 cells from each batch-treatment pairs
select_cell <- NULL
N_sel <- 500
for(r in 1:P){
  
  cell_sub <- sample(which(btp_ind[,3] == r), N_sel, replace = FALSE)
  
  # order the selected cells by cell type
  cell_order <- order(metadata$celltype[cell_sub], cell_sub, decreasing = FALSE)
  cell_sub <- cell_sub[cell_order]
  
  select_cell <- c(select_cell, cell_sub)
}

logy_obs_sub <- log1p(y[,select_cell])
logy_correct_sub <- log1p(x_corrected[, select_cell])

# color bar
color_bar <- cbind(
  color_by_celltype[metadata$celltype[select_cell]],
  color_by_condition[metadata$ condition[select_cell]],
  color_by_batch[metadata$batch[select_cell]])

# color key
break_y <- seq(0,5.5, length.out = 101)
# colorsChoice <- colorRampPalette(c("#1D78B2","#F2F2F2","#C4012D"))
# treatment_key <- colorsChoice(100)

jpeg(paste0("Images/heatmap_logy_obs_v",ver,".jpg"),
     width = 1440, height = 960)
heat_logmu <- heatmap_effects(logy_obs_sub, 
                              color_bar,
                              break_effects = break_y)
dev.off()

jpeg(paste0("Images/colorkey_logy_obs_v",ver,".jpg"),width = 480, height = 240)
draw_keys(logy_obs_sub,
          color_bar,
          break_effects = break_y)
dev.off()

# break_mean <- seq(min(logmu_est_treat), max(logmu_est_treat), length.out = 8)
jpeg(paste0("Images/heatmap_logy_corrected_v",ver,".jpg"),
     width = 1440, height = 960)
heat_logmu <- heatmap_effects(logy_correct_sub, 
                              color_bar,
                              break_effects = break_y)
dev.off()

###############################################
# UMAP for the corrected count data by DIFseq #
###############################################
# Draw PCA and UMAP for corrected read counts
sim_DIFseq_obj <- CreateSeuratObject(counts = x_corrected, meta.data = metadata, project = "sim_DIFseq") 

# scaling
sim_DIFseq_obj <- ScaleData(sim_DIFseq_obj)

# Run PCA
sim_DIFseq_obj <- RunPCA(sim_DIFseq_obj, features = rownames(x_corrected))
ElbowPlot(sim_DIFseq_obj, ndims= 50)

# Run UMAP
sim_DIFseq_obj <- RunUMAP(sim_DIFseq_obj, dims = 1:15)

PCA_DIFseq <- sim_DIFseq_obj[["pca"]]@cell.embeddings
UMAP_DIFseq <- sim_DIFseq_obj[["umap"]]@cell.embeddings

metadata_df$UMAP1 <- UMAP_DIFseq[,1]
metadata_df$UMAP2 <- UMAP_DIFseq[,2]

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

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_shuffled.pdf"),width = 8, height = 8)
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
    # legend.text = element_text(size = 24),
    # legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2))
print(p)
dev.off()

## Color by cell type
color_by <- "CellType"
color_group <- color_by_celltype

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_shuffled.pdf"),width = 8, height = 8)
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
    # legend.text = element_text(size = 24),
    # legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2))
print(p)
dev.off()

## Color by cell type
color_by <- "Condition"
color_group <- color_by_condition

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_shuffled.pdf"),width = 12, height = 8)
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
    # legend.text = element_text(size = 24),
    # legend.title = element_text(size = 32, face="bold"),
    panel.background = element_rect(colour = "black",size = 2))
print(p)
dev.off()


##########################
# Differential abundance #
##########################
### DA across conditions
# Load the MC samples of cell type labels and 
# switch labels to be consistent with true cell type labels
w_MC <- read.table(paste0(dir_inference,"w_MC.txt"))
label_switch <- as.numeric(names(sort(est_switch)))
w_MC <- matrix(label_switch[unlist(w_MC)], ncol = 10)

# calculate the Louis estimator of $\pi_{sk}$ and 
# then conduct differential abundance inference
# function to obtain the Louis estimators of cell type proportion Pi_{sk}
Louis_variance_MCESM <- function(cell_proportion, cell_labels_MC){
  
  if(is.list(cell_proportion)){
    cell_proportion <- unlist(cell_proportion)
  }
  
  n_rep <- ncol(cell_labels_MC)
  K <- length(cell_proportion)
  Km1 <- length(cell_proportion) - 1
  
  sum_w <- matrix(NA, K, n_rep)
  for(r in 1:n_rep){
    for(k in 1:K){
      sum_w[k,r] <- sum(cell_labels_MC[,r]==k)
    }
  }
  
  # Add 0.01 for singular values
  sum_w <- sum_w + 0.01
  
  #print(sum_w)
  
  first_der <- rep(0, Km1)
  first_der_cross <- matrix(0, Km1, Km1)
  second_der <- matrix(0, Km1, Km1)
  
  for(r in 1:n_rep){
    temp_first <- sum_w[1:(K - 1),r]/cell_proportion[1:(K-1)] - sum_w[K]/cell_proportion[K]
    first_der <- first_der + temp_first
    first_der_cross <- first_der_cross + temp_first %*% t(temp_first)
    second_der <- second_der - 
      diag(sum_w[1:(K - 1),r]/(cell_proportion[1:(K-1)])^2) - 
      sum_w[K,r]/(cell_proportion[K])^2
  }
  
  first_der <- first_der/n_rep
  first_der_cross <- first_der_cross/n_rep
  second_der <- second_der/n_rep
  
  obs_variance <- -second_der + first_der_cross - first_der %*% t(first_der)
  obs_variance <- solve(obs_variance)
  
  return(obs_variance)
}

Diff.Abundance <- function(.pi, w, meta, ref = ncol(.pi), dim = c("Cond", "Batch", "Pair"), subset = NULL){
  
  .K <- ncol(.pi)
  .S <- nrow(.pi)
  s_infor <- meta$sample
  .cname <- colnames(.pi)
  .cname <- .cname[-ref]
  # browser()
  
  if(dim == "Cond"){
      # get the indices of each cell and the number of treatments
      d_infor <- factor(meta$condition)
    }else if(dim == "Batch"){
      d_infor <- factor(meta$batch)
    }else if(dim == "Pair"){
      d_infor <- factor(meta$pair)
    }
    
    if(is.null(subset)){

      .dname <- levels(d_infor)
      .d <- length(.dname)
      .dsub <- 1:.d

    }else{

      .dsub <- subset
      .d <- length(subset)
      .dname <- levels(d_infor)[subset]

    }

    # Sample_set stores samples belonging to each treatment
    Sample_set <- list()
    for (j in 1:.d) {
      Sample_set[[j]] <- unique(s_infor[as.numeric(d_infor) == .dsub[j]])
    }

    res <- matrix(NA, .d * (.d - 1) / 2, 2 * .K)
    row_names <- NULL
    
    # compute Louis' estimator
lvar <- array(NA, dim = c(.S, .K - 1, .K - 1))

# put the reference to the last column
if(ref != .K){
  .pi <- cbind(.pi[,-ref] , .pi[,ref])
  pos_ref <- which(w == ref)
  pos_lat <- which(w > ref)
  w[pos_ref] <- .K
  w[pos_lat] <- w[pos_lat] - 1
}

for(s in 1:.S){
  cell_index <- which(s_infor==s)
  lvar[s,,] <- Louis_variance_MCESM(.pi[s,],w[cell_index,]) 
}

index <- 1
for (j1 in 1:(.d - 1)) {
  for (j2 in (j1+1):.d) {
    
    set1 <- Sample_set[[j1]]
    set2 <- Sample_set[[j2]]
    
    ns1 <- length(set1)
    ns2 <- length(set2)
    
    # Difference in proportion
    if(ns1 > 1){
      pi_mean_1 <- apply(.pi[set1, ], 2, mean)
    }else{
      pi_mean_1 <- .pi[set1, ]
    }
    
    if(ns2 > 1){
      pi_mean_2 <- apply(.pi[set2, ], 2, mean)
    }else{
      pi_mean_2 <- .pi[set2, ]
    }
    
    pi_dif <- pi_mean_1 - pi_mean_2
    
    # Compute within-group covariance
    Var_within <- matrix(0, .K - 1, .K - 1)
    for (s in set1) {
      Var_within <- Var_within + lvar[s, ,] / ns1 ^ 2
    }
    
    for (s in set2) {
      Var_within <- Var_within + lvar[s, ,] / ns2 ^ 2
    }
    
    # Compute between-group covariance
    Var_est <- Var_within
    
    if(ns1 > 1){
      Var_est <- Var_est + cov(.pi[set1, 1:(.K - 1)])/ns1
    }
    
    if(ns2 > 1){
      Var_est <- Var_est + cov(.pi[set2, 1:(.K - 1)])/ns2
    }
    
    # test statistics
    stat_celltype <- pi_dif[1:(.K - 1)] / sqrt(diag(Var_est))
    stat_overall <- t(pi_dif[1:(.K - 1)]) %*% solve(Var_est) %*% pi_dif[1:(.K - 1)]
    
    # pval
    p_celltype <- pnorm(abs(stat_celltype), lower.tail = FALSE) * 2
    p_overall <- pchisq(stat_overall, df = .K - 1, lower.tail = FALSE)
    
    # # adjust by BH
    # p_adj <- p.adjust(p_celltype, method = "BH")
    
    res[index, 2 * 1:.K - 1] <- c(stat_celltype, stat_overall)
    res[index, 2 * 1:.K] <-  c(p_celltype, p_overall)    
    row_names <- c(row_names, paste(.dname[j1], "vs", .dname[j2]))
    
    index <- index + 1
  }
}


rownames(res) <- row_names
colnames(res) <- paste(rep(c(.cname,"Overall"),each = 2),rep(c("stat", "p-val"),.K),sep = ":")

return(res)
}

colnames(pi_ordered) <- paste("Type",1:K_opt)

# Batch 1, 2, 3 in Condition 1
DA_Batch_Cond1 <- Diff.Abundance(pi_ordered, w_MC, meta = metadata, 
                                 ref = 5, dim = "Pair", subset = c(1,3,5,7))
xtable(DA_Batch_Cond1, digits = 3)

# Batch 1 in treatment 1, 2 and 3
DA_Cond_Batch1 <- Diff.Abundance(pi_ordered, w_MC, meta = metadata, 
                                 ref = 5, dim = "Pair", subset = c(1,2))
# round(DA_batch1_treat123, digits = 4)
xtable(DA_Cond_Batch1, digits = 3)

# Across three treatments over all the batches
DA_Cond <- Diff.Abundance(pi_ordered, w_MC, meta = metadata, 
                          ref = 5, dim = "Cond", subset = c(1,2))
xtable(DA_Cond, digits = 3)

color_by <- "Stats"
DA_Cond_ref5 <- Diff.Abundance(pi_ordered, w_MC, meta = metadata, 
                               ref = 5, dim = "Cond", subset = c(1,2))

DA_Cond_ref1 <- Diff.Abundance(pi_ordered, w_MC, meta = metadata, 
                               ref = 1, dim = "Cond", subset = c(1,2))

Stat_cond <- matrix(NA, Num_Cond, K_opt)
Stat_cond[1,] <- c(DA_Cond_ref5[c(1,3,5,7)], DA_Cond_ref1[7])
Stat_cond[2,] <- -Stat_cond[1,]

# Build signed -log10(pval)
sign_cond <- sign(Stat_cond)
logpval_cond <- - sign_cond * (log10(2) + 
                                 pnorm(abs(Stat_cond), lower.tail = FALSE, log.p = TRUE)/log(10))


metadata_df$log10_pval <- NULL
metadata_df$DA_Stats <- NULL
for(t in 1:Num_Cond){
  for(k in 1:K_opt){
    index_interest <- (metadata$celltype == k) & (metadata$condition == t)
    metadata_df$DA_Stats[index_interest] <- Stat_cond[t,k]
    metadata_df$log10_pval[index_interest] <- logpval_cond[t,k]
  }
}

type <- "UMAP"
feature1 <- "UMAP1"
feature2 <- "UMAP2"
color_by <- "log10_pval"

metadata_df$log10_pval <- ifelse(metadata_df$log10_pval > 20, 20, metadata_df$log10_pval)
metadata_df$log10_pval <- ifelse(metadata_df$log10_pval < -20, -20, metadata_df$log10_pval)

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,".pdf"),
    width = 8, height = 8)
p <- ggplot(metadata_df, aes_string(x = feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_gradient2(low = "#2c7fb8", mid = "#e0e0e0",high = "#f03b20", midpoint = 0) +
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
  labs(color = "Log10(pval)")
print(p)
dev.off()

pdf(paste0("Images/",type,"_",proj,"_",method,"_v",ver,"_by_",color_by,"_with_legend.pdf"),
    width = 8, height = 8)
p <- ggplot(metadata_df, aes_string(x = feature1, y= feature2, color = color_by)) +
  geom_point(size=0.5) + theme_classic() +
  scale_color_gradient2(low = "#2c7fb8", mid = "#e0e0e0",high = "#f03b20", midpoint = 0) +
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
  labs(color = "Log10(pval)")
print(p)
dev.off()

############################
# Differential expression #
###########################
.fdrDEindicator <- function(xi, kappa){
  
  ind_intr <- xi <= kappa
  fdr <- sum(xi[ind_intr])/sum(ind_intr)
  
  return(fdr)
}

# Calculate the DE posterior probability threshold
.postprob_DE_thr_fun <- function(xi, fdr_threshold=0.05){
  
  kappa_fdr_matr <- NULL
  kappa_set <- sort(unique(xi))
  
  kappa_ind <- which(kappa_set < 0.5 & kappa_set > fdr_threshold)
  
  for(i in kappa_ind){
    
    kappa <- kappa_set[i]
    fdr <- .fdrDEindicator(xi, kappa=kappa)
    
    if(fdr > fdr_threshold){
      break
    }
  }
  
  kappa <- kappa_set[i-1]
  return(kappa)
}

# Estimate intrinsic gene indicators
.estimate_IG_indicators <- function(xi, postprob_DE_threshold = 0.5){
  
  EstL <- xi
  EstL[xi >= postprob_DE_threshold] <- 0
  EstL[xi <= postprob_DE_threshold] <- 1
  # message("The output format is a matrix.\n")
  # message(paste0("Each row represents a gene, and each column",
  #               " corresponds to a cell type from 2 to K\n"))
  return(EstL)
}

# Intrinsic gene index
.IG_index <- function(EstIGindicators){
  ind <- which(rowSums(EstIGindicators) > 0)
  message(c(length(ind), " intrinsic genes are found.\n"))
  message("The output format is a vector implying the intrinsic gene",
          " indices.\n")
  return(ind)
}

# Load log(Pr(L = 0)) and log(Pr(J = 0))
PrL_est <- matrix(unlist(read.table(paste0(dir_inference,"PrL_est.txt"))), nrow = G)
PrJ_est <- matrix(unlist(read.table(paste0(dir_inference,"PrJ_est.txt"))), nrow = G)

PrL_est <- exp(PrL_est)
PrJ_est <- exp(PrJ_est)

colnames(PrL_est) <- paste0("K",1:K_opt)
colnames(PrJ_est) <- paste0("C",rep(1:Num_Cond, each = K_opt),"K",rep(1:K_opt, Num_Cond))

pval_col <- cbind(PrL_est[,-1], PrJ_est[, K_opt + 1:((Num_Cond - 1) * K_opt)])

xi_thres <- .postprob_DE_thr_fun(pval_col, 0.05)
DE_est <- .estimate_IG_indicators(pval_col, xi_thres)

D_est <- apply(DE_est[,1:(K_opt-1)], 1, sum) > 0
num_intri <- sum(D_est)

E_est <- apply(DE_est[,K_opt - 1 + 1:((Num_Cond - 1) * K_opt)], 1, sum) > 0
num_DE <- sum(E_est)

### Intrinsic genes across cell types

# `PrL_est.txt` stores the log posterior probability of $\beta_{gk}$ 
# not far away from zero, that is, 
# $\log(\text{Pr}(L_{gk}=0))$, and `D_est` stores the identified intrinsic genes 
# by controlling FDR at level 0.05. As we know the true intrinsic genes, 
# we can compute the true FDR and FNR.  
logPrL_est <- read.table(paste0(dir_inference,"PrL_est.txt"))
logPrL_est <- matrix(unlist(logPrL_est),G,K_opt)

D_syn <- apply(beta.syn != 0, 1, sum) > 0
D_est <- unlist(read.table(paste0(dir_inference,"D_est.txt")))

## FDR
message("The false discovery rate of intrinsic genes is ",sum(D_est == 1 & D_syn == 0)/sum(D_est == 1))

## FNR
message("The false negative rate of intrinsic genes is ",sum(D_est == 0 & D_syn == 1)/sum(D_est == 0))

################################################
# consider the expression proportions of genes #
################################################
thres_adj <- 0.3
## intrinsic genes
L_est <- DE_est[,1:(K_opt - 1)]
L_adj <- matrix(0, G, K_opt - 1)

expr_prop <- array(NA, dim = c(G, K_opt, 2))

for(k in 2:K_opt){
  celltype_k <- which(w_est == k)
  for(g in which(L_est[,k-1] == 1)){
    
    yg <- y[g, ]
    
    ygk <- yg[celltype_k]
    ygnk <- yg[-celltype_k]
    
    expr_prop[g,k,1] <- sum(ygk > 0)/length(ygk)
    expr_prop[g,k,2] <- sum(ygnk > 0)/length(ygnk)
    
    if(expr_prop[g,k,1] > thres_adj & expr_prop[g,k,2] > thres_adj){
      L_adj[g, k - 1] <- 1 
    }
  }
}

D_adj <- apply(L_adj, 1, sum) > 0
sum(D_adj)

# By controlling FDR at 0.05, 
# we perfectly identify all intrinsic genes. 

### Cell-type-specific DE genes
# Similarly, `PrJ_est.txt` stores the log posterior probability of $\eta_{tgk}$ 
# not far away from zero, that is, $\log(\text{Pr}(J_{tgk}=0))$. 

# est_switch <- auto_switch(w, w_est, k)
logPrJ_est <- read.table(paste0(dir_inference,"PrJ_est.txt"))
logPrJ_est <- matrix(unlist(logPrJ_est),G,K_opt * Num_Cond)

eta_ordered <- matrix(NA, G,K_opt * Num_Cond)
logPrJ_ordered <- matrix(NA, G,K_opt * Num_Cond)
for(t in 1:Num_Cond){
  for(k in 1:K_opt){
    logPrJ_ordered[,(t-1) * K_opt + k] <- logPrJ_est[,(t-1) * K_opt + est_switch[k]]
    eta_ordered[,(t-1) * K_opt + k] <- eta_est[,(t-1) * K_opt + est_switch[k]]
  }
}


control_genes <- unlist(read.table(paste0(dir_inference,"control_genes.txt")))
all(eta.syn[control_genes==1,] == 0)

E_syn <- matrix(0,G,Num_Cond * K_opt)
E_syn <- eta.syn != 0

E_est <- matrix(unlist(read.table(paste0(dir_inference,"E_est.txt"))),G,K_opt * Num_Cond)
E_ordered <- E_est
for(t in 1:Num_Cond){
  E_ordered[, (t - 1) * K_opt + 1:K_opt] <- E_est[,(t - 1) * K_opt + est_switch] == 1
}

# Exclude the reference condition and control genes and only keep valid DE tests
col_valid <- 1:((Num_Cond - 1) * K_opt) + K_opt

E_syn_valid <- E_syn[control_genes==0, col_valid]
E_est_valid <- E_ordered[control_genes==0, col_valid]

## FDR
message("The false discovery rate of cell-type-specific DE genes is ",sum(E_est_valid == 1 & E_syn_valid == 0)/sum(E_est_valid == 1))

## FNR
message("The false negative rate of cell-type-specific DE genes is ",sum(E_est_valid == 0 & E_syn_valid == 1)/sum(E_est_valid == 0))

J_est <- DE_est[,K_opt - 1 + 1:((Num_Cond - 1) * K_opt)]
J_adj <- matrix(0, G, (Num_Cond - 1) * K_opt)

expr_prop_condition <- array(NA, dim = c(G, Num_Cond,K_opt, 2))

for(t in 2:Num_Cond){
  for(k in 1:K_opt){
    celltype_t <- which(w_est == k & metadata$condition == t)
    celltype_other <- which(w_est == k & metadata$condition != t)
    for(g in which(J_est[,k + (t-2) * K_opt] == 1)){
      
      yg <- y[g, ]
      
      ygt <- yg[celltype_t]
      ygnt <- yg[celltype_other]
      
      expr_prop_condition[g,t,k,1] <- sum(ygt > 0)/length(ygt)
      expr_prop_condition[g,t,k,2] <- sum(ygnt > 0)/length(ygnt)
      print(expr_prop_condition[g,t,k,])
      
      if(expr_prop_condition[g,t,k,1] > thres_adj & expr_prop_condition[g,t,k,2] > thres_adj){
        J_adj[g, k + (t-2) * K_opt] <- 1 
      }
    }
  }
}

# Number of DE gene for each cell type
apply(J_est, 2, sum)
apply(J_adj, 2, sum)

E_adj <- apply(J_adj, 1, sum) > 0
sum(E_adj)

# Draw ROC curve for cell-type-specific DE genes
logPrJ_est_valid <- logPrJ_ordered[control_genes==0, col_valid]
logPrJ_est_valid <- matrix(logPrJ_est_valid, nrow = (N - sum(control_genes)) * K_opt)
logPrJ_est_valid <- apply(logPrJ_est_valid, 1, min)

DE_ROC_DIFseq <- data.frame(DE = as.vector(ifelse(E_syn_valid, "DE", "Not DE")),
                            logPrJ = logPrJ_est_valid)

DE_ROC_DIFseq$logPrJ <- ifelse(DE_ROC_DIFseq$logPrJ < -100, -100, DE_ROC_DIFseq$logPrJ)
write.table(DE_ROC_DIFseq, paste0("Output/",method,"_cell_type_specific_DE_genes_for_ROC_v",ver,".txt"))
