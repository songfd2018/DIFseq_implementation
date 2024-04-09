##########################################
# Analyze the Covid-19 dataset by DIFseq #
##########################################
rm(list=ls())

library(pheatmap)
library(ggplot2)
library(viridis)
library(dplyr)

proj <- "covid"
method <- "DIFseq"
ver <- 1
dir <- "Covid/"
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

adt_count <- read.table(paste0("RawCountData/adt_data_",proj,"_v",ver,".txt"))
adt_temp <- as.matrix(adt_count[,-1])
rownames(adt_temp) <- adt_count[,1]
adt_count <- adt_temp

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
k_sel <- 8:15
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
# load cell type labels
w_est <- read.table(paste0(dir_inference,"w_est.txt"))
w_est <- unlist(w_est)
# table(w_est, metadata$CellType)

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

# Load the imputed read counts
x_imputed <- read.table(paste0(dir_inference,"imputed_count.txt"))
colnames(x_imputed) <- paste0(metadata$Sample_id, metadata$Barcode)
rownames(x_imputed) <- gene_list

#############################
# Annotate cell type labels #
#############################
est_switch <- c(12, 9, 3, 11, 4, 5, 6, 8, 2, 10, 1, 7)
celltype_name <-  c("B cells", "Plasma cells", "CD4+ T cells 1", 
                    "CD4+ T cells 2", "CD8+ T cells", 
                    "NK 1", "NK 2", "pDC", "cDC", 
                    "CD14+ Mono",  "CD16+ Mono", "Platelet")

pi_ordered <- pi_est[,est_switch]
colnames(pi_ordered) <- celltype_name
rownames(pi_ordered) <- metadata$Severity[!duplicated(metadata$Sample_id)]

##########################################
# Draw Heatmap of ADT levels by pheatmap #
##########################################
num_adt <- nrow(adt_count)

adt_celltype <- matrix(NA, num_adt, K_opt)

for(k in 1:K_opt){
  adt_celltype[,k] <- apply(log1p(adt_count[,w_est == est_switch[k]]), 1, mean)
}
range(adt_celltype)

colnames(adt_celltype) <- colnames(pi_ordered)
rownames(adt_celltype) <- rownames(adt_count)

pdf(paste0("Images/heatmap_scaled_LogAvgAdt_Celltype_v",ver,"_K",K_opt,".pdf"), width = 6, height = 8)
pheatmap(adt_celltype, 
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         scale = "row") 
dev.off()

##############################################
# Draw cell type proportions for each sample #
##############################################
colnames(pi_ordered) <- celltype_name

sample_id_ordered <- c("259","279","280","258","265",
                       "nCOV1EUHM","nCOV7EUHM","nCOV0029EUHM",
                       "nCOV3EUHM", "nCOV6EUHM","nCOV021EUHM","nCOV024EUHM")

pi.ordered.df <- data.frame(
  Sample=factor(rep(unique(metadata$Sample_id),K_opt),
                levels=sample_id_ordered),
  CellType = factor(rep(colnames(pi_ordered), each = S), levels = celltype_name),
  Prop = as.vector(pi_ordered))

color_by_celltype <- c("#fb8072", # For B cells
                       "#EA9A7E", # For Plasma cells
                       "#fdb462", # For CD4+ T cells 1
                       "#F4CB56", # For CD4+ T cells 2
                       "#ffed6f", # For CD8+ T cells
                       "#b3de69", # For NK 1 cells
                       "#69ED9B", # For NK 2 cells
                       "#4AE2CC", # For pDC
                       "#4BDDE1", # For cDC
                       "#bebada", # For CD14+ Monocyte
                       "#cab2d6", # For CD16+ Monocyte
                       "#bc80bd"# For platelet
)

pdf(file = paste0("Images/CellTypeProportion_",proj,"_v",ver,".pdf"), width = 14, height = 8)
p <- ggplot(data = pi.ordered.df, mapping = aes(x = Sample, fill = CellType, y = Prop)) + 
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

### output intrinsic genes for KEGG pathway analysis
Intrinsic_genes <- gene_list[D_adj==1]
write.table(Intrinsic_genes, file = paste0("Output/Intrinsic_gene_",proj,"_v",ver,".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

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

# Differentially expressed genes across conditions
E_combinded <- J_adj[,1:K_opt] | J_adj[,1:K_opt + K_opt]
E_ordered <- E_combinded[, est_switch]
colnames(E_ordered) <- colnames(pi_ordered)
intri_acorss_cond <- apply(E_ordered,2,sum)


#gene_list <- unlist(gene_list)
for(k in 1:K_opt){
  Intrinsic_genes <- gene_list[E_ordered[,k]==1]
  write.table(Intrinsic_genes, file = paste0("Output/DE_gene_",proj,"_v",ver,"_",celltype_name[k],".txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

####################################################
# Barplot of number of cell-type-specific DE genes #
####################################################
J_ordered <- J_adj[,c(est_switch, est_switch+K_opt)]
color_by_condition <- c("#ccebc5","#b3cde3","#fbb4ae")

df_num_DE <- data.frame(CellType = rep(celltype_name,2),
                        Condition = rep(c("Moderate", "Severe"),each = K_opt),
                        Num_DE_gene = apply(J_ordered, 2, sum))


pdf(file = paste0("Images/DE_gene_number_",proj,"_v",ver,".pdf"), width = 8, height = 9)
p <- ggplot(data = df_num_DE, mapping = aes(x = CellType, fill = Condition, y = log(Num_DE_gene))) + 
  geom_col(color = "#939393", position=position_dodge(), width =0.75) +
  scale_fill_manual(values=color_by_condition[2:3]) +
  scale_x_discrete(limits=rev) +
  coord_flip() +
  labs(x = "Cell Type", y = "Log(No. of DE genes)",fill = "Condition") +
  theme(axis.text=element_text(size=20), 
        axis.title=element_text(size=24,face="bold"),
        # axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.title = element_text(size=20, face = "bold"), #change legend title font size
        legend.text = element_text(size=20))
p
dev.off()

######################################################
# Draw the heatmap of condition effects by heatmap.3 #
######################################################
## Draw the heatmap of condition effects by `heatmap3.R` and `pheatmap`
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

condition_effects_by_celltype <- eta_organized[E_adj,]

# range
break_combined <- seq(-10,10, length.out = 101)
colorsChoice <- colorRampPalette(c("#1D78B2","#F2F2F2","#C4012D"))
color_key <- colorsChoice(100)

# Color bar
color_bar <- cbind(rep(color_by_condition,K_opt),
                   rep(color_by_celltype, each = Num_Cond))

jpeg(paste0("Images/heatmap_condition_effects_by_celltype_v",ver,"_K",K_opt,".jpg"), width = 480, height = 960, quality = 100)
heat_logmu <- heatmap_effects(condition_effects_by_celltype,
                              color_bar,
                              break_effects = break_combined,
                              labRow = FALSE, labCol = FALSE,
                              color_key = color_key)
dev.off()

jpeg(paste0("Images/colorkey_condition_effects_by_celltype_v",ver,".jpg"),width = 480, height = 240)
draw_keys(condition_effects_by_celltype, 
          color_key = color_key,
          break_effects = break_combined)
dev.off()

###################################
# Differential abundance analysis #
###################################
w_MC <- read.table(paste0(dir_inference,"w_MC.txt"))
names(est_switch) <- 1:K_opt
label_switch <- as.numeric(names(sort(est_switch)))
w_MC <- matrix(label_switch[unlist(w_MC)], ncol = 10)
pi_ordered <- pi_est[,est_switch]


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

colnames(pi_ordered) <- celltype_name
meta_da <- data.frame(batch = metadata$Batch, 
                      condition = factor(metadata$Severity), 
                      pair = btp_ind[,3], 
                      sample = btp_ind[,4])


# Two batches across condition 1, 2 and 3
DA_cond <- Diff.Abundance(pi_ordered, w_MC, meta = meta_da, ref = K_opt, dim = "Cond")

####################################
# Obtain the corrected read counts #
####################################
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

start_time<-Sys.time()
message("Calculate corrected read counts:")
x_corrected <-adjusted_values(y_obs, btp_ind, K_opt,
                              alpha_est, beta_est, eta_est,nu_est,delta_est,phi_est,w_est)
write.table(x_corrected, file = paste0("Inference_K",K_opt,"/corrected_count.txt"), row.names = FALSE, col.names = FALSE)
end_time<-Sys.time()
running_time<-difftime(end_time, start_time, units = "mins")
message("It takes ",running_time," mins to calculate the corrected data.")

rownames(metadata) <- paste0(metadata$Sample_id, metadata$Barcode)

##################################################
## Dot plot of the marker gene expression levels #
##################################################
marker_genes <- c("IFI27", "IFITM3", "ISG15")

# Adjust cell type labels
names(est_switch) <- 1:K_opt
label_switch <- as.numeric(names(sort(est_switch)))
w_ordered <- label_switch[w_est]

df_marker <- NULL

for(gene in marker_genes){
  
  log_mean_counts <- rep(NA, K_opt * Num_Cond)
  expr_prop <- rep(NA, K_opt * Num_Cond)
  
  corrected_count <- x_corrected[gene_list == gene, ]
  
  df_count <- data.frame(DIFseq_type = celltype_name[w_ordered],
                         Condition = metadata$Severity,
                         Count = corrected_count)
  
  
  df_summary <- df_count %>% 
    group_by(DIFseq_type, Condition) %>%
    summarise(LogMeanExp = mean(log1p(Count)), 
              ExprProp = mean(Count > 0))
  
  
  df_summary$Marker <- gene
  
  df_marker <- rbind(df_marker, df_summary)
}

df_marker$DIFseq_type <- factor(df_marker$DIFseq_type, levels = celltype_name)

pdf(file = paste0("Images/Dotplot_IFN_genes_",proj,"_v",ver,".pdf"), width = 10, height = 10)
ggplot(df_marker, aes(x = DIFseq_type, y = Condition)) + 
  geom_point(aes(color = LogMeanExp, size = ExprProp), alpha = 0.8) + facet_grid(Marker ~.) + 
  # scale_y_discrete(limits=rev) +
  scale_color_gradient(low = "#ffeda0", high = "#e31a1c") + 
  labs(x = "Cell Type", y = "Condition",color = "Average Expression", size = "Precent Expressed") +
  theme(axis.text=element_text(size=20), 
        axis.title=element_text(size=24,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.title = element_text(size=20, face = "bold"), 
        legend.text = element_text(size=20),
  ) + 
  guides(fill=guide_legend(title="-log10(pval)"))
dev.off()
