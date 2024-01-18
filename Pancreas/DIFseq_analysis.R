rm(list=ls())
# Analyze pancreas datasets by DIFseq
library(mclust)
library(ggplot2)
library(viridis)
library(pheatmap)
library(RColorBrewer)
library(Seurat)
library(dplyr)
library(patchwork)

proj <- "pancreas"
method <- "DIFseq"
ver <- 1
dir <- "Pancreas/"
setwd(dir)
source("heatmap3.R")

if(!dir.exists("Images")){
  dir.create("Images")
}

if(!dir.exists("Output")){
  dir.create("Output")
}

####################################
# Load raw count data and metadata #
####################################
y_obs <- read.table(paste0("RawCountData/count_data_",proj,"_v",ver,".txt"))
y_obs <- as.matrix(y_obs)

# Load dimension
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

# Load metadata
metadata <- read.table(paste0("RawCountData/metadata_",proj,"_v",ver,".txt"), header = TRUE)
# colnames(metadata) <- c("Batch", "Condition", "CellType")

# Load the gene list
gene_list <- unlist(read.table(paste0("RawCountData/gene_list_",proj,"_v",ver,".txt"),stringsAsFactors = F))

# Add row and column names for the raw count data matrix
rownames(y_obs) <- gene_list
colnames(y_obs) <- metadata$Sample

# indices of cells
btp_ind <- read.table(paste0("RawCountData/tbinfor_",proj,"_v",ver,".txt"))

##################################################
# Select the optimal number of cell types by BIC #
##################################################
k_sel <- 3:12
BIC.record <- rep(NA, length(k_sel))
for (k in k_sel) {
  # load BIC
  inference_dir <- paste0("Inference_K",k,"/")
  BIC <- unlist(read.table(paste0(inference_dir, "BIC.txt")))
  BIC.record[k - k_sel[1] + 1] <- BIC
}

# Draw the BIC plot
pdf(paste0("Images/BIC_plot_",proj,"_v",ver,".pdf"),width = 6, height = 8)
par(mar = c(5.1,6.1,4.1,2.1)) 
plot(k_sel,BIC.record,xlab= "K",ylab = "BIC",type="n",cex.axis=2,cex.lab=3)
points(k_sel,BIC.record,type="b",pch=19,cex=3)
dev.off()

# Consider the optimal number of cell types
K_opt <- k_sel[which.min(BIC.record)]
dir_inference <- paste0("Inference_K",K_opt,"/")

#######################################
# Load the estimated parameter values #
#######################################
# logistic parameters 
gamma_est <- read.table(paste0(dir_inference,"gamma_est.txt"))

# load cell type labels
w_est <- read.table(paste0(dir_inference,"w_est.txt"))
w_est <- unlist(w_est)
table(w_est, metadata$CellType)

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
tau0 <- 0.1
tau1_est <- ptau1[1]
pbeta_est <- ptau1[2]
peta_est <- ptau1[3]

#############################
# Annotate cell type labels #
#############################
annotated_names <- c("Alpha", "Beta", "Delta", "Gamma", "Acinar", "Ductal", "other")
metadata$CellType <- factor(metadata$CellType, levels = annotated_names)
table_celltype <- table(w_est, metadata$CellType)
est_switch <- c(1,2,6,5,4,3)
table_celltype <- table_celltype[est_switch,]

celltype_name <- c("Alpha", "Beta", "Delta&Gamma", "Acinar", "Ductal", "Others")
DIFseq_labels <- w_est
for(i in 1:K_opt){
  DIFseq_labels[DIFseq_labels == est_switch[i]] <- celltype_name[i]
}
DIFseq_labels <- factor(DIFseq_labels, levels = celltype_name)

ARI_DIFseq <- adjustedRandIndex(w_est, metadata$CellType)
write.table(ARI_DIFseq, file = paste0("Output/ARI_",method,"_",proj,"_v",ver,".txt"), quote = FALSE)

#############################################
# Draw cell type proprotions of each sample #
#############################################
pi_ordered <- pi_est[,est_switch]
colnames(pi_ordered) <- celltype_name
sample_id <- unique(metadata$Donor)
rownames(pi_ordered) <- sample_id

sample_id <- unique(metadata$Donor)

df_prop <- data.frame(Sample=factor(rep(sample_id, K_opt),levels = sample_id),
                      CellType = factor(rep(celltype_name, each = S), levels = celltype_name),
                      Prop = as.vector(pi_ordered))

color_by_condition <- c("#ffffb3","#fbb4ae")
color_by_batch <- viridis(4)
color_by_celltype <- c( "#fb8072","#fdb462","#C7EF85","#8dd3c7","#80b1d3","#bebada")

pdf(file = paste0("Images/CellTypeProportion_",proj,"_v",ver,".pdf"), width = 12, height = 8)
p <- ggplot(data = df_prop, mapping = aes(x = Sample, fill = CellType, y = Prop)) + 
  geom_col(width = 1, color = "#939393") +
  scale_fill_manual(values=color_by_celltype) +
  labs(x = "Sample", y = "Cell Type Proportion",colour = "Cell Type") +
  theme(axis.text=element_text(size=20), 
        axis.title=element_text(size=24,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.title = element_text(size=20, face = "bold"), #change legend title font size
        legend.text = element_text(size=20),
        panel.grid = element_blank(),
        panel.background =element_blank(),
        panel.border = element_blank()) +
  guides(fill=guide_legend(title="Cell Type")) +
  ylab("Cell-type proportions")
p
dev.off()

#########################################
# Draw color annotation bar by pheatmap #
#########################################
# Define the column names and row names
pi_pheat <- t(pi_ordered)
colnames(pi_pheat) <- sample_id
annotation_samples <- data.frame(
  Condition = factor(metadata$Disease[!duplicated(metadata$Donor)]),
  Batch = factor(metadata$Study[!duplicated(metadata$Donor)]))
rownames(annotation_samples) <- sample_id

# Specify the color of the annotation bar


names(color_by_condition) <- c("ND", "T2D")
names(color_by_batch) <- unique(metadata$Study)

sample_colors = list(Condition = color_by_condition,
                     Batch = color_by_batch)

pdf(paste0("Images/heatmap_celltype_proportion_v",ver,"_K",K_opt,".pdf"), width = 12, height = 8)
pheatmap(pi_pheat, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdYlBu")))(100),
         breaks = seq(0,1,length.out = 100),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_col = annotation_samples,
         annotation_colors = sample_colors,
         show_colnames = FALSE,
         scale = "none") 
dev.off()

###################################
# Differential abundance analysis #
###################################
w_MC <- read.table(paste0(dir_inference,"w_MC.txt"))
names(est_switch) <- 1:K_opt
label_switch <- as.numeric(names(sort(est_switch)))
w_MC <- matrix(label_switch[unlist(w_MC)], ncol = 10)

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
    # get the indices of each cell and the number of conditions
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
  
  # Sample_set stores samples belonging to each condition
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
      p_celltype <- 2 * (1 - pnorm(abs(stat_celltype)))
      p_overall <- pchisq(stat_overall, df = .K - 1, lower.tail = FALSE)
      
      # # adjust by BH
      # p_adj <- p.adjust(p_celltype, method = "BH")
      
      res[index, 2 * 1:.K - 1] <- c(stat_celltype, stat_overall)
      # res[index, 2 * 1:.K] <-  c(p_adj, p_overall)
      res[index, 2 * 1:.K] <-  c(p_celltype, p_overall)
      
      row_names <- c(row_names, paste(.dname[j1], "vs", .dname[j2]))
      
      index <- index + 1
    }
  }
  
  
  rownames(res) <- row_names
  colnames(res) <- paste(rep(c(.cname,"Overall"),each = 2),rep(c("stat", "pval"),.K),sep = ":")
  
  return(res)
}

meta_da <- data.frame(batch = btp_ind[,1], condition = btp_ind[,2], pair = btp_ind[,3], sample = btp_ind[,4])

# Diff.Abundance(pi_ordered, w_MC, meta = metadata, ref = 1, dim = "Cond", subset = c(1,2), alg = "MCESM")

# Batch 3 in condition 1 and 2
DA_GSE86473 <- Diff.Abundance(pi_ordered[,1:3], w_MC, meta = meta_da, ref = 3, dim = "Pair", subset = c(3,4))
DA_GSE86473

# Exclude the 40th donor HP1526901T2D 
meta_da[meta_da[,4] == 40,3] <- 100

# Batch 4 in condition 1, 2
DA_EMTAB5061 <- Diff.Abundance(pi_ordered, w_MC, meta = meta_da, ref = K_opt, dim = "Pair", subset = c(5,6))
DA_EMTAB5061

####################################
# Differential expression analysis #
####################################
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
colnames(PrJ_est) <- paste0("T",rep(1:Num_Cond, each = K_opt),"K",rep(1:K_opt, Num_Cond))

pval_col <- cbind(PrL_est[,-1], PrJ_est[, K_opt + 1:((Num_Cond - 1) * K_opt)])

xi_thres <- .postprob_DE_thr_fun(pval_col, 0.01)
DE_est <- .estimate_IG_indicators(pval_col, xi_thres)

D_est <- apply(DE_est[,1:(K_opt-1)], 1, sum) > 0
num_intri <- sum(D_est)

E_est <- apply(DE_est[,K_opt - 1 + 1:((Num_Cond - 1) * K_opt)], 1, sum) > 0
num_DE <- sum(E_est)

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
    
    yg <- y_obs[g, ]
    
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

### output intrinsic genes for pathway analysis
Intrinsic_genes <- gene_list[D_adj==1]
write.table(Intrinsic_genes, file = paste0("Intrinsic_gene_",proj,"_v",ver,".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

J_est <- DE_est[,K_opt - 1 + 1:((Num_Cond - 1) * K_opt)]
J_adj <- matrix(0, G, (Num_Cond - 1) * K_opt)

expr_prop_condition <- array(NA, dim = c(G, Num_Cond,K_opt, 2))

for(t in 2:Num_Cond){
  for(k in 1:K_opt){
    celltype_t <- which(w_est == k & btp_ind[,2] == t)
    celltype_other <- which(w_est == k & btp_ind[,2] != t)
    for(g in which(J_est[,k + (t-2) * K_opt] == 1)){
      
      yg <- y_obs[g, ]
      
      ygt <- yg[celltype_t]
      ygnt <- yg[celltype_other]
      
      expr_prop_condition[g,t,k,1] <- sum(ygt > 0)/length(ygt)
      expr_prop_condition[g,t,k,2] <- sum(ygnt > 0)/length(ygnt)
      # print(expr_prop_condition[g,t,k,])
      
      if(expr_prop_condition[g,t,k,1] > thres_adj & expr_prop_condition[g,t,k,2] > thres_adj){
        J_adj[g, k + (t-2) * K_opt] <- 1 
      }
    }
  }
}

E_adj <- apply(J_adj, 1, sum) > 0
sum(E_adj)

# # Differentially expressed genes across conditions
E_combinded <- J_adj[,1:K_opt]
E_ordered <- E_combinded[, est_switch]
colnames(E_ordered) <- celltype_name
intri_acorss_cond <- apply(E_ordered,2,sum)

celltype_names <- colnames(pi_ordered)

#gene_list <- unlist(gene_list)
for(k in 1:K_opt){
  Intrinsic_genes <- gene_list[E_ordered[,k]==1]
  write.table(Intrinsic_genes, file = paste0("Output/DE_gene_",proj,"_v",ver,"_",celltype_names[k],".txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

#########################################
# Draw the heatmap of condition effects #
#########################################
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

heatmap_effects <- function(effects_matrix, color_bar, color_key = NULL, break_effects = NULL, ...){
  
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
                       dendrogram = "row",#with cluster tree
                       Rowv = TRUE, Colv = FALSE,
                       # labRow = FALSE, labCol = FALSE,
                       ColSideColors = color_bar,
                       lmat=rbind(c(5,4),c(0,1), c(3,2)),#1=heatmap, 2=row dendogram, 3=col dendogram, 4= key
                       lhei=c(0.3,0.4,3.6),
                       lwid = c(0.3,3),
                       col=color_key, breaks = break_effects, key = FALSE, ...)
  
  return(out_fig)
  
}

eta_ordered <- eta_est
for(t in 2:Num_Cond){
  eta_ordered[,(t-1) * K_opt + 1:K_opt] <- eta_est[, (t-1) * K_opt + est_switch]
}

eta_organized <- eta_ordered
for(k in 1:K_opt){
  eta_organized[,(k-1) * Num_Cond + 1:Num_Cond] <- eta_ordered[,0:(Num_Cond-1) * K_opt + k]
}

rownames(eta_organized) <- gene_list
condition_effects_by_celltype <- eta_organized[E_adj,]

####################
# Draw by pheatmap #
####################
# Color bar
cond_id <- paste(
  rep(celltype_name, each = Num_Cond), 
  rep(1:Num_Cond, K_opt), 
  sep = "_")

colnames(condition_effects_by_celltype) <- cond_id

annotation_cond <- data.frame(
  Condition = factor(rep(c("ND","T2D"), K_opt)),
  CellType = factor(rep(celltype_name, each = Num_Cond), 
                    levels = celltype_name))

rownames(annotation_cond) <- cond_id

# Specify the color of the annotation bar
# color_by_condition <- c("#ccebc5","#b3cde3","#fbb4ae")
# color_by_celltype <- rainbow(K_opt)

names(color_by_celltype) <- celltype_name

cond_colors = list(Condition = color_by_condition,
                   CellType = color_by_celltype)

# logmu_ordered <- alpha_est + beta_ordered + eta_ordered
pdf(paste0("Images/pheatmap_condition_effects_by_celltype_v",ver,"_K",K_opt,".pdf"), width = 8, height = 12)
pheatmap(condition_effects_by_celltype,
         color = colorRampPalette(c("#1D78B2","#F2F2F2","#C4012D"))(101),
         breaks = seq(-4,4,0.08),
         cluster_cols = FALSE,
         annotation_col = annotation_cond,
         annotation_colors = cond_colors,
         show_colnames = FALSE,
         scale = "none") 
dev.off()

#####################
# Draw by heatmap.3 #
#####################
# range
break_combined <- seq(-3,3, length.out = 101)
colorsChoice <- colorRampPalette(c("#1D78B2","#F2F2F2","#C4012D"))
color_key <- colorsChoice(100)

# Color bar
color_bar <- cbind(rep(color_by_condition,K_opt),
                   rep(color_by_celltype, each = Num_Cond))

pdf(paste0("Images/heatmap_condition_effects_by_celltype_v",ver,"_K",K_opt,".pdf"), width = 6, height = 6)
heat_logmu <- heatmap_effects(condition_effects_by_celltype,
                              color_bar,
                              break_effects = break_combined,
                              labCol = FALSE,
                              color_key = color_key)
dev.off()

jpeg(paste0("Images/colorkey_condition_effects_by_celltype_v",ver,".jpg"),width = 480, height = 240)
draw_keys(condition_effects_by_celltype, 
          color_key = color_key,
          break_effects = break_combined)
dev.off()

#########################################################################
# Correct batch effects and visualize the corrected read counts by UMAP #
#########################################################################
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
    
    if(i %% 100 == 0){
      print(paste("Finish the correction of", i, "cells..."))
      
    }
  }
  return(CorrectedCount)
}

x_imputed <- read.table(paste0(dir_inference,"imputed_count.txt"))
colnames(x_imputed) <- metadata$Sample
rownames(x_imputed) <- gene_list

start_time<-Sys.time()
message("Calculate corrected read counts:")
x_corrected <-adjusted_values(x_imputed, btp_ind, K_opt,
                              alpha_est, beta_est, eta_est,nu_est,delta_est,phi_est,w_est)
end_time<-Sys.time()
running_time<-difftime(end_time, start_time, units = "mins")
message("It takes ",running_time," mins to calculate the corrected data.")

rownames(metadata) <- metadata$Sample

###############################################
# UMAP for the corrected count data by DIFseq #
###############################################
# Draw PCA and UMAP for corrected read counts
pancreas_DIFseq_obj <- CreateSeuratObject(counts = x_corrected, 
                                          meta.data = metadata, 
                                          project = "pancreas_DIFseq") 

# scaling
pancreas_DIFseq_obj <- ScaleData(pancreas_DIFseq_obj)

# Run PCA
pancreas_DIFseq_obj <- RunPCA(pancreas_DIFseq_obj, 
                              features = rownames(x_corrected))
ElbowPlot(pancreas_DIFseq_obj, ndims= 50)

# Run UMAP
pancreas_DIFseq_obj <- RunUMAP(pancreas_DIFseq_obj, dims = 1:20)

PCA_DIFseq <- pancreas_DIFseq_obj[["pca"]]@cell.embeddings
UMAP_DIFseq <- pancreas_DIFseq_obj[["umap"]]@cell.embeddings

metadata$UMAP1 <- UMAP_DIFseq[,1]
metadata$UMAP2 <- UMAP_DIFseq[,2]

# Shuffle cells
set.seed(123)
cell_shuffle <- sample(1:N, N)
metadata_shuffled <- metadata[cell_shuffle, ]

# Draw UMAP plots
type <- "UMAP"
feature1 <- "UMAP1"
feature2 <- "UMAP2"

## Color by batch
color_by <- "Study"
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
n.celltype <- length(unique(metadata$CellType))
color_by_celltype_FACS <- brewer.pal(n.celltype, "Set3")
color_group <- color_by_celltype_FCAS

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
color_by <- "Disease"
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

############################
# Dot plot of marker genes #
############################
marker_genes <- c("GCG", "INS", "PPY", "SST", "PRSS1", "KRT19")

# Adjust cell type labels
w_ordered <- label_switch[w_est]

# gene <- "ISG15"

df_marker <- NULL
for(gene in marker_genes){
  log_mean_counts <- rep(NA, K_opt * Num_Cond)
  expr_prop <- rep(NA, K_opt * Num_Cond)
  
  corrected_count <- unlist(x_corrected[gene_list == gene, ])
  
  df_count <- data.frame(DIFseq_type = celltype_name[w_ordered],
                         Condition = metadata$Disease,
                         Count = corrected_count)
  
  df_summary <- df_count %>% 
    group_by(DIFseq_type, Condition) %>%
    summarise(LogMeanExp = mean(log1p(Count)), 
              ExprProp = mean(Count > 0))
  
  
  df_summary$Marker <- gene
  
  df_marker <- rbind(df_marker, df_summary)
}

df_marker$DIFseq_type <- factor(df_marker$DIFseq_type, levels = celltype_name)
df_marker$Marker <- factor(df_marker$Marker, levels = marker_genes)

pdf(file = paste0("Images/Dotplot_marker_genes_",proj,"_v",ver,".pdf"), width = 8, height = 12)
ggplot(df_marker, aes(x = DIFseq_type, y = Condition)) + 
  geom_point(aes(color = LogMeanExp, size = ExprProp), alpha = 0.8) + facet_grid(Marker ~.) + 
  # scale_y_discrete(limits=rev) +
  scale_color_gradient(low = "#F2F2F2", high = "#C4012D") + 
  labs(x = "Cell Type", y = "Condition",color = "Average Expression", size = "Precent Expressed") +
  theme(axis.text=element_text(size=20), 
        axis.title=element_text(size=24,face="bold"),
        strip.text.y = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.title = element_text(size=20, face = "bold"), 
        legend.text = element_text(size=20),
  )
dev.off()

######################################################################
# Draw the violin plot of cell-type-specific DE genes for Beta cells #
######################################################################
# celltype_name <- c("Alpha", "Beta", "Delta&Gamma", "Acinar", "Ductal", "Others")
k <- 2
Intrinsic_beta <- gene_list[E_ordered[,k]==1]
plots_DIFseqMarker <- VlnPlot(pancreas_DIFseq_obj, 
                              features = Intrinsic_beta, 
                              split.by = "Disease", group.by = "CellType",
                              pt.size = 0, combine = FALSE)
wrap_plots(plots = plots_DIFseqMarker, ncol = 2)
