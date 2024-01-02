# Program: 01ModelSelection_pancreas_v53.R
# Select the optimal number of cell types by BIC
# Contents: 
library(ggplot2)

rm(list=ls())
proj <- "covid_sc"
# K <- 5
ver <- 53

# Working directory
# setwd(paste0("D:/CUHKSZ/BUTseq/code/",proj,"/v",ver,"/"))
setwd(paste0("E:/Study/DIFseq/",proj,"/v",ver,"/"))

#################################################
# Select the optimal number of cell type by BIC #
#################################################
load(paste0("Model_selection_post_v",ver,".RData"))

# BIC
write.table(summary_results[[1]], "clipboard", sep="\t", row.names=FALSE)

# ARI_1
write.table(summary_results[[2]], "clipboard", sep="\t", row.names=FALSE)

# ARI_2
write.table(summary_results[[3]], "clipboard", sep="\t", row.names=FALSE)

# Iteration number
write.table(summary_results[[4]], "clipboard", sep="\t", row.names=FALSE)

# Running time
write.table(summary_results[[5]]/3600, "clipboard", sep="\t", row.names=FALSE)

# The number of intrinsic genes across cell types
write.table(summary_results[[6]], "clipboard", sep="\t", row.names=FALSE)

# The number of cell-type-specific DE genes
write.table(summary_results[[7]], "clipboard", sep="\t", row.names=FALSE)

# Estimated slab variance
write.table(summary_results[[8]], "clipboard", sep="\t", row.names=FALSE)

# Draw BIC plot
df_BIC <- data.frame(CellTypeNum = 3:9,
                     BIC = apply(summary_results[[1]][-7,],1,min,na.rm = TRUE))

pdf(paste0("Images/BIC_plot_pancreas_v",ver,".pdf"),width = 8, height = 6)
ggplot(data=df_BIC, aes(x=CellTypeNum, y=BIC)) +
  geom_line()+
  geom_point()
dev.off()

# ###################################
# # Compare metadata of v36 and v37 #
# ###################################
# metadata_v36 <- read.table("../v36/RawCountData/metadata_covid_sc_v36.txt", header = TRUE)
# metadata_v37 <- read.table("../v37/RawCountData/metadata_covid_sc_v37.txt", header = TRUE)
#   
# mismatch <- which(metadata_v36$Barcode != metadata_v37$Barcode)


################
# Trend of ARI #
K_opt <- 6
r_opt <- 3
dir_inference <- paste0("K",K_opt,"_r",r_opt,"/Inference_K",K_opt,"/")
time_con <- read.table(file = paste0(dir_inference,"time_consumption.txt"),header = FALSE)
colnames(time_con) <- c("E step for spike-and-slab priors","M step for pi",
                        "MCE step","SM step for gamma","SM step for abnep",
                        "M step for delta","M step for p and tau",
                        "Calculate the likelihood")
time_per_iter <- apply(time_con,2 ,mean)
time_percentage <- round(time_per_iter/sum(time_per_iter) * 100, digits = 2)
time_percentage

iter_infor <- unlist(read.table(paste0(dir_inference,"iter_infor.txt")))

iter_num <- nrow(time_con)
plot(c(1,iter_num), c(0,500),type = "n")
lines(1:iter_num, time_con[,3], col = 3)
lines(1:iter_num, time_con[,4], col = 4)
lines(1:iter_num, time_con[,5], col = 5)
abline(v = iter_infor[3], lty = 2, col = 2)
abline(v = sum(iter_infor[3:4]), lty = 2, col = 2)
abline(v = sum(iter_infor[3:5]), lty = 2, col = 2)

## Take a look at the convergence
Traceplot <- function(par, subset = NULL, range = NULL, 
                      ...){
  N_iter <- nrow(par)
  Dim <- ncol(par)
  if(!is.null(subset)){
    par <- par[,subset]
  }
  
  # determine the range of the plot
  if(is.null(range)){
    for(i in 1:N_iter){
      range_i <- range(par[i,])
      range <- range(c(range,range_i))
    }
  }
  
  p <- plot(c(1,N_iter), range, type = "n", 
            xlab = "IterNum", ...)
  for(g in 1:Dim){
    p <- lines(par[,g],col = g)
  }
  return(p)
}

iter_dir <- paste0("K",K_opt,"_r",r_opt,"/MCESM_iter_K",K_opt,"/")

mini_batch <- 1000
sub_trace <- paste("Traceplot of minibatch size equal to",mini_batch) 

alpha_post <- read.big.matrix(paste0(iter_dir,"alpha_iter.txt"), type = "double", sep = " ")
jpeg(filename = paste0("Images/Traceplot_alpha_",proj,"_v",ver,"_K",K_opt,"_r",r_opt,".jpg"), width = 800, height = 600)
Traceplot(alpha_post, ylab = expression(alpha[g]), 
          main = sub_trace)
dev.off()

beta_post <- read.big.matrix(paste0(iter_dir,"beta_iter.txt"), type = "double", sep = " ")
jpeg(filename = paste0("Images/Traceplot_beta_",proj,"_v",ver,"_K",K_opt,"_r",r_opt,".jpg"), width = 800, height = 600)
Traceplot(beta_post, ylab = expression(beta[gk]), 
          main = sub_trace)
dev.off()

eta_post <- read.big.matrix(paste0(iter_dir,"eta_iter.txt"), type = "double", sep = " ")
jpeg(filename = paste0("Images/Traceplot_eta_",proj,"_v",ver,"_K",K_opt,"_r",r_opt,".jpg"), width = 800, height = 600)
Traceplot(eta_post, ylab = expression(eta[tgk]), 
          main = sub_trace)
dev.off()

nu_post <- read.big.matrix(paste0(iter_dir,"nu_iter.txt"), type = "double", sep = " ")
jpeg(filename = paste0("Images/Traceplot_nu_",proj,"_v",ver,"_K",K_opt,"_r",r_opt,".jpg"), width = 800, height = 600)
Traceplot(nu_post, ylab = expression(nu[bg]),
          main = sub_trace)
dev.off()

delta_post <- read.big.matrix(paste0(iter_dir,"delta_iter.txt"), type = "double", sep = " ")
jpeg(filename = paste0("Images/Traceplot_delta_",proj,"_v",ver,"_K",K_opt,"_r",r_opt,".jpg"), width = 800, height = 600)
Traceplot(delta_post, ylab = expression(delta[btgi]), 
          main = sub_trace)
dev.off()

phi_post <- read.big.matrix(paste0(iter_dir,"phi_iter.txt"), type = "double", sep = " ")
jpeg(filename = paste0("Images/Traceplot_phi_",proj,"_v",ver,"_K",K_opt,"_r",r_opt,".jpg"), width = 800, height = 600)
Traceplot(phi_post, ylab = expression(phi[bg]), 
          main = sub_trace)
dev.off()

# ################
# # trace of ARI #
# ################
# library(bigmemory)
# library(mclust)
# 
# # K6_r11
# K <- 7
# rep <- 1
# w_iter <- read.big.matrix(paste0("K",K,"_r",rep,"/MCESM_iter_K",K,"/w_iter.txt"), sep = " ", type = "integer")
# Record_num <- nrow(w_iter)
# ARI_iter <- rep(NA, Record_num)
# 
# # load metadata
# metadata <- read.table(paste0("../v49/RawCountData/metadata_pancreas_v49.txt"))
# w_true <- metadata[,3]
# 
# for(i in 1:Record_num){
#   ARI_iter[i] <- adjustedRandIndex(w_true,w_iter[i,])
# }
# 
# iter_infor <- unlist(read.table(paste0("K",K,"_r",rep,"/Inference_K",K,"/iter_infor.txt")))
# iter_step1 <- iter_infor[3] + 1
# iter_step2 <- sum(iter_infor[3:4]) + 1
# 
# pdf(paste0("Images/ARI_trace_",proj,"_v",ver,"_K",K,"_r",rep,".pdf"), width = 8, height = 6)
# plot(1:Record_num * 5, ARI_iter, type = "l",
#      xlab = "IterNum", ylab = "ARI",
#      main = "ARI trace of the pancreas study")
# abline(v = iter_step1, lty = 2, col = 2)
# abline(v = iter_step2, lty = 2, col = 2)
# dev.off()
# 
# w_prev <- unlist(read.table("../v49/K6_r5/Inference_K6/w_est.txt"))
# w_cur <- unlist(read.table("K6_r12/Inference_K6/w_est.txt"))
# w_cur2 <- unlist(read.table("K6_r11/Inference_K6/w_est.txt"))
# 
# adjustedRandIndex(w_prev,w_cur)
# adjustedRandIndex(w_prev,w_true)
# adjustedRandIndex(w_true,w_cur)
# 
# table(w_prev, w_cur)
# 
# table(w_cur, w_true)
# 
# table(w_cur2, w_true)
# 
# table(w_prev, w_true)
# 
# # check the likelihood calculation
# rep <- 2
# K <- 6
# dim <-unlist(read.table(paste0("RawCountData/dim_pancreas_v",ver,".txt")))
# N <- dim[1]
# G <- dim[2]
# B <- dim[3]
# NumTreat <- dim[4]
# P <- dim[5]
# 
# gamma.est <- unlist(read.table(paste0("K",K,"_r",rep,"/Inference_K",K,"/gamma_est.txt")))
# gamma.est <- matrix(gamma.est, B, 2)
# 
# pi.est <- unlist(read.table(paste0("K",K,"_r",rep,"/Inference_K",K,"/pi_est.txt")))
# pi.est <- matrix(pi.est, P, K)
# 
# alpha.est <- unlist(read.table(paste0("K",K,"_r",rep,"/Inference_K",K,"/alpha_est.txt")))
# 
# beta.est <- unlist(read.table(paste0("K",K,"_r",rep,"/Inference_K",K,"/beta_est.txt")))
# beta.est <- matrix(beta.est, G, K)
# 
# eta.est <- unlist(read.table(paste0("K",K,"_r",rep,"/Inference_K",K,"/eta_est.txt")))
# eta.est <- matrix(eta.est, G, K * NumTreat)
# 
# nu.est <- unlist(read.table(paste0("K",K,"_r",rep,"/Inference_K",K,"/nu_est.txt")))
# nu.est <- matrix(nu.est, G, B)
# 
# delta.est <- unlist(read.table(paste0("K",K,"_r",rep,"/Inference_K",K,"/delta_est.txt")))
# 
# phi.est <- unlist(read.table(paste0("K",K,"_r",rep,"/Inference_K",K,"/phi_est.txt")))
# phi.est <- matrix(nu.est, G, B)
# 
# y <- unlist(read.table(paste0("RawCountData/count_data_pancreas_v",ver,".txt")))
# y <- matrix(y, G, N)
# 
# bt_infor <- read.table(paste0("RawCountData/tbinfor_pancreas_v",ver,".txt"))
# 
# # calculate observed likelihood
# 
# for(i in 1:N){
#   
#   b <- bt_infor[i, 1]
#   t <- bt_infor[i, 2]
#   p <- bt_infor[i, 3]
#   
#   # calculate sum_g Pr(y_btig|\Theta) for each k
#   log_ratio <- rep(0, K)
#   for(g in 1:G){
#     log_mu <- alpha.est[g] + beta.est[g,] + eta.est[g, (t-1) * K + 1:K] + nu.est[g,b] + delta.est[i]
#     p_temp <- exp(log_mu)/(exp(log_mu) + phi.est[g,b])
#     
#   }
#   
# }
# 
# PrL.est <- read.table("../v70/K6_r6/Inference_K6/PrL_est.txt")
# head(PrL.est)
# 
# vec_PrL <- unlist(PrL.est[,2:6])
# thres <- sort(vec_PrL)
# 
# for(i in 1:length(thres)){
#   FDR <- sum(vec_PrL[vec_PrL < thres[i]])/* 
# }
# 
# PrL.est <- read.table("../v73/K6_r6/Inference_K6/PrL_est.txt")

###############
# Convergence #
###############
