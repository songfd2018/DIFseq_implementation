# Program: 01ModelSelection_pancreas_v53.R
# Select the optimal number of cell types by BIC
# Contents: 
library(ggplot2)

rm(list=ls())
proj <- "pancreas"
# K <- 5
ver <- 161

# Working directory
# setwd(paste0("D:/CUHKSZ/BUTseq/code/",proj,"/v",ver,"/"))
setwd(paste0("E:/Study/DIFseq/",proj,"/v",ver,"/"))

#################################################
# Select the optimal number of cell type by BIC #
#################################################
load(paste0("Model_selection_post_v",ver,".RData"))

# BIC
write.table(summary_results[[1]], "clipboard", sep="\t", row.names=FALSE)

# ARI
write.table(summary_results[[2]], "clipboard", sep="\t", row.names=FALSE)

# Iteration number
write.table(summary_results[[3]], "clipboard", sep="\t", row.names=FALSE)

# Running time
write.table(summary_results[[4]]/3600, "clipboard", sep="\t", row.names=FALSE)

# The number of intrinsic genes across cell types
write.table(summary_results[[5]], "clipboard", sep="\t", row.names=FALSE)

# The number of cell-type-specific DE genes
write.table(summary_results[[6]], "clipboard", sep="\t", row.names=FALSE)

# Estimated slab variance
write.table(summary_results[[7]], "clipboard", sep="\t", row.names=FALSE)

# Draw BIC plot

df_BIC <- data.frame(CellTypeNum = 3:9,
                     BIC = apply(summary_results[[1]][-7,],1,min,na.rm = TRUE))

pdf(paste0("Images/BIC_plot_pancreas_v",ver,".pdf"),width = 8, height = 6)
ggplot(data=df_BIC, aes(x=CellTypeNum, y=BIC)) +
  geom_line()+
  geom_point()
dev.off()

########################################################################
# Compare different types of BIC and quantify the approximation of BIC #
########################################################################
# Load observed values
dim <-unlist(read.table(paste0("RawCountData/dim_",proj,"_v",ver,".txt")))
N <- dim[1]
S <- dim[2]
G <- dim[3]
B <- dim[4]
NumTreat <- dim[5]
P <- dim[6]

y <- as.matrix(read.table(paste0("RawCountData/count_data_",proj,"_v",ver,".txt")))

btp <- read.table(paste0("RawCountData/tbinfor_",proj,"_v",ver,".txt"))
b_ind <- btp[,1]
t_ind <- btp[,2]
s_ind <- btp[,4]

# Load the estimated parameter values
rep <- 5
K <- 6

gamma.est <- unlist(read.table(paste0("K",K,"_r",rep,"/Inference_K",K,"/gamma_est.txt")))
gamma.est <- matrix(gamma.est, B, 2)

pi.est <- unlist(read.table(paste0("K",K,"_r",rep,"/Inference_K",K,"/pi_est.txt")))
pi.est <- matrix(pi.est, S, K)

alpha.est <- unlist(read.table(paste0("K",K,"_r",rep,"/Inference_K",K,"/alpha_est.txt")))

beta.est <- unlist(read.table(paste0("K",K,"_r",rep,"/Inference_K",K,"/beta_est.txt")))
beta.est <- matrix(beta.est, G, K)

eta.est <- unlist(read.table(paste0("K",K,"_r",rep,"/Inference_K",K,"/eta_est.txt")))
eta.est <- matrix(eta.est, G, K * NumTreat)

nu.est <- unlist(read.table(paste0("K",K,"_r",rep,"/Inference_K",K,"/nu_est.txt")))
nu.est <- matrix(nu.est, G, B)

delta.est <- unlist(read.table(paste0("K",K,"_r",rep,"/Inference_K",K,"/delta_est.txt")))

phi.est <- unlist(read.table(paste0("K",K,"_r",rep,"/Inference_K",K,"/phi_est.txt")))
phi.est <- matrix(phi.est, G, B)

w.est <- unlist(read.table(paste0("K",K,"_r",rep,"/Inference_K",K,"/w_est.txt")))
x_imputed <- as.matrix(read.table(paste0("K",K,"_r",rep,"/Inference_K",K,"/imputed_count.txt")))

# Focus on Cell 1
i <- 1
b <- b_ind[i]
t <- t_ind[i]
y_cur <- y[,i]
x_cur <- x_imputed[,i]

logmu <- alpha.est + beta.est + eta.est[, (t-1) * K + 1:K] + nu.est[,b] + delta.est[i]

## Complete data log-likelihood
### Dropout logistic terms
complete_loglike <- sum((gamma.est[b,1] + gamma.est[b,2] * x_cur) * (y_cur > 0) - (log(1 + exp(gamma.est[b,1] + gamma.est[b,2] * x_cur))) * (x_cur > 0))

### NB terms
k_cur <- w.est[i]
complete_loglike <- complete_loglike + sum(lgamma(x_cur + phi.est[,b]) - lgamma(x_cur + 1) - lgamma(phi.est[,b])) 
complete_loglike <- complete_loglike + sum(x_cur * logmu[,k_cur] + phi.est[,b] * log(phi.est[,b]) - (x_cur + phi.est[,b]) * log(phi.est[,b] + exp(logmu[,k_cur])))


#############################
## Observed log-likelihood ##
#############################
upper_ratio <- seq(2,4, by = 0.1)
length_ratio <- length(upper_ratio)
obs_loglike_ratio <- rep(0, length_ratio)

for(j in 1:length_ratio){
  
  # obs_loglike <- 0
  # for(i in 1:N){
  for(i in 1:100){
    
    b <- b_ind[i]
    t <- t_ind[i]
    s <- s_ind[i]
    y_cur <- y[,i]
    
    logmu <- alpha.est + beta.est + eta.est[, (t-1) * K + 1:K] + nu.est[,b] + delta.est[i]
    
    
    loglike_k <- rep(0, K)
    
    for(g in 1:G){
      if(y_cur[g] > 0){
        ## Logistic model
        loglike_k <- loglike_k - log(1 + exp(gamma.est[b,1] + gamma.est[b,2] * y_cur[g]))
        
        ## NB parts
        loglike_k <- loglike_k + lgamma(y_cur[g] + phi.est[g,b]) - lgamma(y_cur[g] + 1) - lgamma(phi.est[g,b])
        loglike_k <- loglike_k + y_cur[g] * logmu[g,] + phi.est[g,b] * log(phi.est[g,b]) - (phi.est[g,b] + y_cur[g]) * log(phi.est[g,b] + exp(logmu[g,]))
        
      }else{
        log_cum <- rep(0, K)
        
        for(k in 1:K){
          
          x_max <- round(upper_ratio[j] * exp(logmu[g,k] + gamma.est[b,1]))
          log_temp <- rep(0, x_max + 1) 
          
          ### For x = 0
          if(x_max > 0){
            log_temp[1] <- phi.est[g,b] * log(phi.est[g,b] / (phi.est[g,b] + exp(logmu[g,k])))
            log_temp[2:(x_max+1)] <-  - log(1+exp( - gamma.est[b,1] - gamma.est[b,2] * (1:x_max))) # logistic model
            log_temp[2:(x_max+1)] <-  log_temp[2:(x_max+1)] + lgamma(1:x_max + phi.est[g,b]) - lgamma(1:x_max + 1) + lgamma(phi.est[g,b]) # normalizing constant
            log_temp[2:(x_max+1)] <-  log_temp[2:(x_max+1)] + logmu[g,k] * (1:x_max) + phi.est[g,b] * log(phi.est[g,b]) - (1:x_max + phi.est[g,b]) * log(exp(logmu[g,k]) + phi.est[g,b]) # 
          }else{
            log_temp <- phi.est[g,b] * log(phi.est[g,b] / (phi.est[g,b] + exp(logmu[g,k])))
          }
          
          log_cum[k] <- max(log_temp) + log(sum(exp(log_temp - max(log_temp))))
        }
        
        loglike_k <- loglike_k + log_cum
        # # Draw the traceplot
        # plot(c(0, x_max), c(0.5,1), type = "n", xlab = "x_max", ylab = "Approximated probability") 
        # 
        # for(k in 1:K){
        #   lines(0:x_max, log_cum[,k], col = k + 1)
        # }
      }
      
      if(sum(is.infinite(loglike_k)) > 0){
        print(g)
        print(loglike_k)
        break
      }
    }
    
    loglike_k <- loglike_k + log(pi.est[s,])
    maxlog_k <- max(loglike_k)
    obs_loglike_ratio[j] <- obs_loglike_ratio[j] + maxlog_k + log(sum(exp(loglike_k - maxlog_k)))
    
    if(i %% 100 == 0){
      print(paste("Finish calculating the log-likelihood of the first",i,"cells."))
    }
  }
}
plot(upper_ratio, obs_loglike_ratio, type = "b")


upper_bound <- 5:10
length_bound <- length(upper_bound)
obs_loglike_bound <- rep(0, length_bound)
for(j in 1:length_bound){
  
  # obs_loglike <- 0
  # for(i in 1:N){
  for(i in 1:100){
    
    b <- b_ind[i]
    t <- t_ind[i]
    s <- s_ind[i]
    y_cur <- y[,i]
    
    logmu <- alpha.est + beta.est + eta.est[, (t-1) * K + 1:K] + nu.est[,b] + delta.est[i]
    
    
    loglike_k <- rep(0, K)
    
    for(g in 1:G){
      if(y_cur[g] > 0){
        ## Logistic model
        loglike_k <- loglike_k - log(1 + exp(gamma.est[b,1] + gamma.est[b,2] * y_cur[g]))
        
        ## NB parts
        loglike_k <- loglike_k + lgamma(y_cur[g] + phi.est[g,b]) - lgamma(y_cur[g] + 1) - lgamma(phi.est[g,b])
        loglike_k <- loglike_k + y_cur[g] * logmu[g,] + phi.est[g,b] * log(phi.est[g,b]) - (phi.est[g,b] + y_cur[g]) * log(phi.est[g,b] + exp(logmu[g,]))
        
      }else{
        log_cum <- rep(0, K)
        
        for(k in 1:K){
          
          x_max <- upper_bound[j]
          log_temp <- rep(0, x_max + 1) 
          
          ### For x = 0
          if(x_max > 0){
            log_temp[1] <- phi.est[g,b] * log(phi.est[g,b] / (phi.est[g,b] + exp(logmu[g,k])))
            log_temp[2:(x_max+1)] <-  - log(1+exp( - gamma.est[b,1] - gamma.est[b,2] * (1:x_max))) # logistic model
            log_temp[2:(x_max+1)] <-  log_temp[2:(x_max+1)] + lgamma(1:x_max + phi.est[g,b]) - lgamma(1:x_max + 1) + lgamma(phi.est[g,b]) # normalizing constant
            log_temp[2:(x_max+1)] <-  log_temp[2:(x_max+1)] + logmu[g,k] * (1:x_max) + phi.est[g,b] * log(phi.est[g,b]) - (1:x_max + phi.est[g,b]) * log(exp(logmu[g,k]) + phi.est[g,b]) # 
          }else{
            log_temp <- phi.est[g,b] * log(phi.est[g,b] / (phi.est[g,b] + exp(logmu[g,k])))
          }
          
          log_cum[k] <- max(log_temp) + log(sum(exp(log_temp - max(log_temp))))
        }
        
        loglike_k <- loglike_k + log_cum
      }
      
      if(max(loglike_k) > 0){
        print(g)
        print(loglike_k)
        break
      }
      
      if(sum(is.infinite(loglike_k)) > 0){
        print(g)
        print(loglike_k)
        break
      }
    }
    
    loglike_k <- loglike_k + log(pi.est[s,])
    maxlog_k <- max(loglike_k)
    obs_loglike_bound[j] <- obs_loglike_bound[j] + maxlog_k + log(sum(exp(loglike_k - maxlog_k)))
    
    if(i %% 100 == 0){
      print(paste("Finish calculating the log-likelihood of the first",i,"cells."))
    }
  }
}
# Figure out the relationship between the fold on mu and the log-likelihood
plot(upper_bound, obs_loglike_bound, type = "b")

i <- 1
b <- b_ind[i]
t <- t_ind[i]
s <- s_ind[i]
y_cur <- y[,i]
w_cur <- w.est[i]
logmu <- alpha.est + beta.est + eta.est[, (t-1) * K + 1:K] + nu.est[,b] + delta.est[i]


# scatter plot of log1p(y_{sig}) and logmu
plot(log1p(y_cur), logmu[,w_cur], xlab = "log(1+y)", ylab = "log(mu hat)")
abline(a=0,b=1, col=2,lty = 2)

# Traceplot of the fold and approximated observed log-likelihood
upper_ratio <- seq(2,12, by = 2)
length_ratio <- length(upper_ratio)
loglike_i_ratio <- rep(0, length_ratio)

for(j in 1:length_ratio){
  loglike_k <- rep(0, K)
  
  for(g in 1:G){
    if(y_cur[g] > 0){
      
      ## Logistic model
      loglike_k <- loglike_k - log(1 + exp(gamma.est[b,1] + gamma.est[b,2] * y_cur[g]))
      
      ## NB parts
      loglike_k <- loglike_k + lgamma(y_cur[g] + phi.est[g,b]) - lgamma(y_cur[g] + 1) - lgamma(phi.est[g,b])
      loglike_k <- loglike_k + y_cur[g] * logmu[g,] + phi.est[g,b] * log(phi.est[g,b]) - (phi.est[g,b] + y_cur[g]) * log(phi.est[g,b] + exp(logmu[g,]))
      
    }else{
      log_cum <- rep(0, K)
      
      for(k in 1:K){
        
        x_max <- round(upper_ratio[j] * exp(logmu[g,k] + gamma.est[b,1]))
        log_temp <- rep(0, x_max + 1) 
        
        ### For x = 0
        if(x_max > 0){
          log_temp[1] <- phi.est[g,b] * log(phi.est[g,b] / (phi.est[g,b] + exp(logmu[g,k])))
          log_temp[2:(x_max+1)] <-  - log(1+exp( - gamma.est[b,1] - gamma.est[b,2] * (1:x_max))) # logistic model
          log_temp[2:(x_max+1)] <-  log_temp[2:(x_max+1)] + lgamma(1:x_max + phi.est[g,b]) - lgamma(1:x_max + 1) - lgamma(phi.est[g,b]) # normalizing constant
          log_temp[2:(x_max+1)] <-  log_temp[2:(x_max+1)] + logmu[g,k] * (1:x_max) + phi.est[g,b] * log(phi.est[g,b]) - (1:x_max + phi.est[g,b]) * log(exp(logmu[g,k]) + phi.est[g,b]) # 
        }else{
          log_temp <- phi.est[g,b] * log(phi.est[g,b] / (phi.est[g,b] + exp(logmu[g,k])))
        }
        
        log_cum[k] <- max(log_temp) + log(sum(exp(log_temp - max(log_temp))))
      }
      
      loglike_k <- loglike_k + log_cum
    }
    
  }
  
  loglike_k <- loglike_k + log(pi.est[s,])
  maxlog_k <- max(loglike_k)
  loglike_i_ratio[j] <- loglike_i_ratio[j] + maxlog_k + log(sum(exp(loglike_k - maxlog_k)))
}

plot(upper_ratio, loglike_i_ratio, type = "b", xlab = "Fold", ylab = paste0("log-likelihood for cell ",i))

# Traceplot of the fixed upper bound and approximated observed log-likelihood
upper_bound <- 5:10
length_bound <- length(upper_bound)
loglike_i_bound <- rep(0, length_bound)

for(j in 1:length_bound){
  loglike_k <- rep(0, K)
  
  for(g in 1:G){
    if(y_cur[g] > 0){
      
      ## Logistic model
      loglike_k <- loglike_k - log(1 + exp(gamma.est[b,1] + gamma.est[b,2] * y_cur[g]))
      
      ## NB parts
      loglike_k <- loglike_k + lgamma(y_cur[g] + phi.est[g,b]) - lgamma(y_cur[g] + 1) - lgamma(phi.est[g,b])
      loglike_k <- loglike_k + y_cur[g] * logmu[g,] + phi.est[g,b] * log(phi.est[g,b]) - (phi.est[g,b] + y_cur[g]) * log(phi.est[g,b] + exp(logmu[g,]))
      
    }else{
      log_cum <- rep(0, K)
      
      for(k in 1:K){
        
        x_max <- upper_bound[j]
        log_temp <- rep(0, x_max + 1) 
        
        ### For x = 0
        if(x_max > 0){
          log_temp[1] <- phi.est[g,b] * log(phi.est[g,b] / (phi.est[g,b] + exp(logmu[g,k])))
          log_temp[2:(x_max+1)] <-  - log(1+exp( - gamma.est[b,1] - gamma.est[b,2] * (1:x_max))) # logistic model
          log_temp[2:(x_max+1)] <-  log_temp[2:(x_max+1)] + lgamma(1:x_max + phi.est[g,b]) - lgamma(1:x_max + 1) - lgamma(phi.est[g,b]) # normalizing constant
          log_temp[2:(x_max+1)] <-  log_temp[2:(x_max+1)] + logmu[g,k] * (1:x_max) + phi.est[g,b] * log(phi.est[g,b]) - (1:x_max + phi.est[g,b]) * log(exp(logmu[g,k]) + phi.est[g,b]) # 
        }else{
          log_temp <- phi.est[g,b] * log(phi.est[g,b] / (phi.est[g,b] + exp(logmu[g,k])))
        }
        
        log_cum[k] <- max(log_temp) + log(sum(exp(log_temp - max(log_temp))))
      }
      
      loglike_k <- loglike_k + log_cum
    }
    
  }
  
  loglike_k <- loglike_k + log(pi.est[s,])
  maxlog_k <- max(loglike_k)
  loglike_i_bound[j] <- loglike_i_bound[j] + maxlog_k + log(sum(exp(loglike_k - maxlog_k)))
}

plot(upper_bound, loglike_i_bound, type = "b", xlab = "x_max", ylab = paste0("log-likelihood for cell ",i))

########################################################################################
# What happens to the approximated log-likelihood when genes have large mu but y_g = 0 #
########################################################################################
which(logmu[,w_cur] > 6 & y_cur == 0)

#################################
# For the first cell in batch 4 #
#################################
i <- 4933
b <- b_ind[i]
t <- t_ind[i]
s <- s_ind[i]
y_cur <- y[,i]
w_cur <- w.est[i]
logmu <- alpha.est + beta.est + eta.est[, (t-1) * K + 1:K] + nu.est[,b] + delta.est[i]


# scatter plot of log1p(y_{sig}) and logmu
plot(log1p(y_cur), logmu[,w_cur], xlab = "log(1+y)", ylab = "log(mu hat)")
abline(a=0,b=1, col=2,lty = 2)

# Traceplot of the fold and approximated observed log-likelihood
upper_ratio <- seq(2,12, by = 2)
length_ratio <- length(upper_ratio)
loglike_i_ratio <- rep(0, length_ratio)

for(j in 1:length_ratio){
  loglike_k <- rep(0, K)
  
  for(g in 1:G){
    if(y_cur[g] > 0){
      
      ## Logistic model
      loglike_k <- loglike_k - log(1 + exp(gamma.est[b,1] + gamma.est[b,2] * y_cur[g]))
      
      ## NB parts
      loglike_k <- loglike_k + lgamma(y_cur[g] + phi.est[g,b]) - lgamma(y_cur[g] + 1) - lgamma(phi.est[g,b])
      loglike_k <- loglike_k + y_cur[g] * logmu[g,] + phi.est[g,b] * log(phi.est[g,b]) - (phi.est[g,b] + y_cur[g]) * log(phi.est[g,b] + exp(logmu[g,]))
      
    }else{
      log_cum <- rep(0, K)
      
      for(k in 1:K){
        
        x_max <- round(upper_ratio[j] * exp(logmu[g,k] + gamma.est[b,1]))
        log_temp <- rep(0, x_max + 1) 
        
        ### For x = 0
        if(x_max > 0){
          log_temp[1] <- phi.est[g,b] * log(phi.est[g,b] / (phi.est[g,b] + exp(logmu[g,k])))
          log_temp[2:(x_max+1)] <-  - log(1+exp( - gamma.est[b,1] - gamma.est[b,2] * (1:x_max))) # logistic model
          log_temp[2:(x_max+1)] <-  log_temp[2:(x_max+1)] + lgamma(1:x_max + phi.est[g,b]) - lgamma(1:x_max + 1) - lgamma(phi.est[g,b]) # normalizing constant
          log_temp[2:(x_max+1)] <-  log_temp[2:(x_max+1)] + logmu[g,k] * (1:x_max) + phi.est[g,b] * log(phi.est[g,b]) - (1:x_max + phi.est[g,b]) * log(exp(logmu[g,k]) + phi.est[g,b]) # 
        }else{
          log_temp <- phi.est[g,b] * log(phi.est[g,b] / (phi.est[g,b] + exp(logmu[g,k])))
        }
        
        log_cum[k] <- max(log_temp) + log(sum(exp(log_temp - max(log_temp))))
      }
      
      loglike_k <- loglike_k + log_cum
    }
    
  }
  
  loglike_k <- loglike_k + log(pi.est[s,])
  maxlog_k <- max(loglike_k)
  loglike_i_ratio[j] <- loglike_i_ratio[j] + maxlog_k + log(sum(exp(loglike_k - maxlog_k)))
}

plot(upper_ratio, loglike_i_ratio, type = "b", xlab = "Fold", ylab = paste0("log-likelihood for cell ",i))

# Traceplot of the fixed upper bound and approximated observed log-likelihood
upper_bound <- 5:10
length_bound <- length(upper_bound)
loglike_i_bound <- rep(0, length_bound)

for(j in 1:length_bound){
  loglike_k <- rep(0, K)
  
  for(g in 1:G){
    if(y_cur[g] > 0){
      
      ## Logistic model
      loglike_k <- loglike_k - log(1 + exp(gamma.est[b,1] + gamma.est[b,2] * y_cur[g]))
      
      ## NB parts
      loglike_k <- loglike_k + lgamma(y_cur[g] + phi.est[g,b]) - lgamma(y_cur[g] + 1) - lgamma(phi.est[g,b])
      loglike_k <- loglike_k + y_cur[g] * logmu[g,] + phi.est[g,b] * log(phi.est[g,b]) - (phi.est[g,b] + y_cur[g]) * log(phi.est[g,b] + exp(logmu[g,]))
      
    }else{
      log_cum <- rep(0, K)
      
      # x_max <- 100
      # plot(c(0,x_max),c(-10,5),type = "n", xlab = "Upper bound", ylab = "log-likelihood")
      
      for(k in 1:K){
        
        x_max <- upper_bound[j]
        # x_max <- 100
        
        log_temp <- rep(0, x_max + 1) 
        
        ### For x = 0
        if(x_max > 0){
          log_temp[1] <- phi.est[g,b] * log(phi.est[g,b] / (phi.est[g,b] + exp(logmu[g,k])))
          log_temp[2:(x_max+1)] <-  - log(1+exp( - gamma.est[b,1] - gamma.est[b,2] * (1:x_max))) # logistic model
          log_temp[2:(x_max+1)] <-  log_temp[2:(x_max+1)] + lgamma(1:x_max + phi.est[g,b]) - lgamma(1:x_max + 1) - lgamma(phi.est[g,b]) # normalizing constant
          log_temp[2:(x_max+1)] <-  log_temp[2:(x_max+1)] + logmu[g,k] * (1:x_max) + phi.est[g,b] * log(phi.est[g,b]) - (1:x_max + phi.est[g,b]) * log(exp(logmu[g,k]) + phi.est[g,b]) # 
        }else{
          log_temp <- phi.est[g,b] * log(phi.est[g,b] / (phi.est[g,b] + exp(logmu[g,k])))
        }
        
        log_cum[k] <- max(log_temp) + log(sum(exp(log_temp - max(log_temp))))
        
        # lines(0:x_max, log_temp, col = k + 1)
      }
      
      loglike_k <- loglike_k + log_cum
    }
    
  }
  
  loglike_k <- loglike_k + log(pi.est[s,])
  maxlog_k <- max(loglike_k)
  loglike_i_bound[j] <- loglike_i_bound[j] + maxlog_k + log(sum(exp(loglike_k - maxlog_k)))
}

plot(upper_bound, loglike_i_bound, type = "b", xlab = "x_max", ylab = paste0("log-likelihood for cell ",i))

# How each term in the infinite sum changes with respect to x
k <- 6
g <- which(logmu[,w_cur] > 6 & y_cur == 0)[2]
x_max <- round(20 * exp(logmu[g,k] + gamma.est[b,1]))
# x_max <- 100

log_temp <- rep(0, x_max + 1) 

### For x = 0
if(x_max > 0){
  log_temp[1] <- phi.est[g,b] * log(phi.est[g,b] / (phi.est[g,b] + exp(logmu[g,k])))
  log_temp[2:(x_max+1)] <-  - log(1+exp( - gamma.est[b,1] - gamma.est[b,2] * (1:x_max))) # logistic model
  log_temp[2:(x_max+1)] <-  log_temp[2:(x_max+1)] + lgamma(1:x_max + phi.est[g,b]) - lgamma(1:x_max + 1) - lgamma(phi.est[g,b]) # normalizing constant
  log_temp[2:(x_max+1)] <-  log_temp[2:(x_max+1)] + logmu[g,k] * (1:x_max) + phi.est[g,b] * log(phi.est[g,b]) - (1:x_max + phi.est[g,b]) * log(exp(logmu[g,k]) + phi.est[g,b]) # 
}else{
  log_temp <- phi.est[g,b] * log(phi.est[g,b] / (phi.est[g,b] + exp(logmu[g,k])))
}

plot(0:x_max, log_temp, type ="l", xlab = "x", ylab = "f(u)")
plot(0:x_max, cumsum(exp(log_temp)), type = "l", xlab = "u", ylab = "g(u)", ylim = c(0,0.3))
abline(v = exp(logmu[g,k] + gamma.est[b,1]), col = 2, lty = 2)
abline(v = 2 * exp(logmu[g,k] + gamma.est[b,1]), col = 3, lty = 2)

# approximated upper bound
q_new <- phi.est[g,b]/(exp(logmu[g,k]) * (1-exp(gamma.est[b,2])) + phi.est[g,b])
sum_upper_bound <- exp(gamma.est[b,1]) * q_new^phi.est[g,b]

# Pr(y = 0)
pry0 <- 1/(1+exp(gamma.est[b,1])) * (phi.est[g,b] / (phi.est[g,b] + exp(logmu[g,k])))^phi.est[g,b]

abline(h = sum_upper_bound, col = 4, lty = 3)

plot(1:x_max, - log(1+exp( - gamma.est[b,1] - gamma.est[b,2] * (1:x_max))))


#########################################
# |gamma| is large, mu is close to zero #
#########################################
loglike_sample <- function(obs = 0, upper = 10, logmu, phi, gamma){
  
  # obs gives the observed count, an integer
  # upper represents the upper bound in the infinite sum, an integer or a vector
  # mu denotes the mean expression levels, it can be a vector
  # phi denotes the over-dispersion parameter, a scalar
  # gamma represents the logistic parameter, a two-component vector
  # output: a matrix with upper + 1 rows and K columns
  
  
  K <- length(logmu)
  # browser()
  
  if(obs > 0){
    
    log_out <- rep(0, K)
    ## Logistic model
    log_out <- log_out - log(1 + exp(gamma[1] + gamma[2] * obs))
    
    ## NB parts
    log_out <- log_out + lgamma(obs + phi) - lgamma(obs + 1) - lgamma(phi)
    log_out <- log_out + obs * logmu + phi * log(phi) - 
      (phi + obs) * log(phi + exp(logmu))
    
  }else{
    
    log_out <- matrix(0, upper + 1, K)
    for(k in 1:K){
      
      log_temp <- rep(0, upper + 1) 
      
      ### For x = 0
      if(upper > 0){
        log_temp[1] <- phi * log(phi / (phi + exp(logmu[k])))
        log_temp[2:(upper+1)] <-  gamma[1] + gamma[2] * (1:upper) - log(1+exp(gamma[1] + gamma[2] * (1:upper))) # logistic model
        log_temp[2:(upper+1)] <-  log_temp[2:(upper+1)] + lgamma(1:upper + phi) - 
          lgamma(1:upper + 1) - lgamma(phi) # normalizing constant
        log_temp[2:(upper+1)] <-  log_temp[2:(upper+1)] + logmu[k] * (1:upper) + 
          phi * log(phi) - (1:upper + phi) * log(exp(logmu[k]) + phi) # 
        
        log_out[,k] <- max(log_temp) + log(cumsum(exp(log_temp - max(log_temp))))
      }else{
        log_out[k] <- phi * log(phi / (phi + exp(logmu)))
      }
      
    }
  }
  
  return(log_out)
}

loglike_upper <- function(obs = 0, logmu, phi, gamma, times = 3){
  
  # obs gives the observed count, an integer
  # mu denotes the mean expression levels, it can be a vector
  # phi denotes the over-dispersion parameter, a scalar
  # gamma represents the logistic parameter, a two-component vector
  # output: approximated log-likelihood
  
  K <- length(logmu)
  # browser()
  
  if(obs > 0){
    
    log_out <- rep(0, K)
    ## Logistic model
    log_out <- log_out - log(1 + exp(gamma[1] + gamma[2] * obs))
    
    ## NB parts
    log_out <- log_out + lgamma(obs + phi) - lgamma(obs + 1) - lgamma(phi)
    log_out <- log_out + obs * logmu + phi * log(phi) - 
      (phi + obs) * log(phi + exp(logmu))
    
  }else{
    
    log_out <- rep(0, K)
    log_pr0 <- phi * (log(phi) - log(phi + exp(logmu)))
    log_out <- log_pr0
    t <- 1
    
    while(t <= times){
      
      logq <- log(phi) - log(exp(logmu) * (1-exp(t * gamma[2])) + phi)
      logHt <- t * gamma[1] + phi * logq
      log0t <- log_pr0 + t * gamma[1]
      
      max_log <- max(log_out, logHt, log0t)
      log_out <- max_log + log(exp(log_out - max_log) + 
                                 (-1)^(t+1) * exp(logHt - max_log) +
                                 (-1)^t * exp(log0t - max_log) )
      
      t <- t + 1
    }
  }
  
  return(log_out)
}

####################################
# Case 1: small mu and small gamma #
####################################
logmu_cur <- -3
upper_cur <- 10
phi_cur <- 1
gamma_cur <-  c(-0.5,-2.5)

loglike1 <- loglike_sample(obs = 0, upper = upper_cur, 
                           logmu = logmu_cur, phi = phi_cur, gamma = gamma_cur)

loglike_upper1 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 1)


loglike_upper2 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 3)

loglike_upper3 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 5)
# Draw the traceplot
plot(loglike1,type = "l", ylim = c(-0.05,-0.04),
     xlab = "x_max", ylab = "loglikelihood",
     main = paste0("logmu = ",logmu_cur,
                   ", gamma_0 = ",gamma_cur[1],
                   ", gamma_1 = ",gamma_cur[2]))
abline(h = loglike_upper1, col = 2, lty = 2)
abline(h = loglike_upper2, col = 3, lty = 2)
abline(h = loglike_upper3, col = 4, lty = 2)
legend("topleft",legend = paste0("t = ",c(1,3,5)), lty =2, col = 2:4)

####################################
# Case 2: large mu and small gamma #
####################################
logmu_cur <- 10
upper_cur <- 1000
phi_cur <- 1
gamma_cur <-  c(-0.5,-2.5)

loglike2 <- loglike_sample(obs = 0, upper = upper_cur, 
                           logmu = logmu_cur, phi = phi_cur, gamma = gamma_cur)

loglike_upper1 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 1)


loglike_upper2 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 3)

loglike_upper3 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 5)
# Draw the traceplot
plot(loglike2,type = "l", ylim = c(-10,-9.9),
     xlab = "x_max", ylab = "loglikelihood",
     main = paste0("logmu = ",logmu_cur,
                   ", gamma_0 = ",gamma_cur[1],
                   ", gamma_1 = ",gamma_cur[2]))
abline(h = loglike_upper1, col = 2, lty = 2)
abline(h = loglike_upper2, col = 3, lty = 2)
abline(h = loglike_upper3, col = 4, lty = 2)
legend("topleft",legend = paste0("t = ",c(1,3,5)), lty =2, col = 2:4)

#######################################
# Case 3: moderate mu and small gamma #
#######################################
logmu_cur <- 1
upper_cur <- 100
phi_cur <- 1
gamma_cur <-  c(-0.5,-2.5)

loglike <- loglike_sample(obs = 0, upper = upper_cur, 
                           logmu = logmu_cur, phi = phi_cur, gamma = gamma_cur)

loglike_upper1 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 1)


loglike_upper2 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 3)

loglike_upper3 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 5)
# Draw the traceplot
plot(loglike,type = "l", ylim = c(-1.5,-1),
     xlab = "x_max", ylab = "loglikelihood",
     main = paste0("logmu = ",logmu_cur,
                   ", gamma_0 = ",gamma_cur[1],
                   ", gamma_1 = ",gamma_cur[2]))
abline(h = loglike_upper1, col = 2, lty = 2)
abline(h = loglike_upper2, col = 3, lty = 2)
abline(h = loglike_upper3, col = 4, lty = 2)
legend("topleft",legend = paste0("t = ",c(1,3,5)), lty =2, col = 2:4)

####################################
# Case 4: small mu and large gamma #
####################################
logmu_cur <- -3
upper_cur <- 100
phi_cur <- 1
gamma_cur <-  c(-0.5,-0.002)

loglike <- loglike_sample(obs = 0, upper = upper_cur, 
                          logmu = logmu_cur, phi = phi_cur, gamma = gamma_cur)

loglike_upper1 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 1)


loglike_upper2 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 3)

loglike_upper3 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 5)
# Draw the traceplot
plot(loglike,type = "l", ylim = c(-0.05,0),
     xlab = "x_max", ylab = "loglikelihood",
     main = paste0("logmu = ",logmu_cur,
                   ", gamma_0 = ",gamma_cur[1],
                   ", gamma_1 = ",gamma_cur[2]))
abline(h = loglike_upper1, col = 2, lty = 2)
abline(h = loglike_upper2, col = 3, lty = 2)
abline(h = loglike_upper3, col = 4, lty = 2)
legend("topleft",legend = paste0("t = ",c(1,3,5)), lty =2, col = 2:4)

#######################################
# Case 5: moderate mu and large gamma #
#######################################
logmu_cur <- 1
upper_cur <- 100
phi_cur <- 1
gamma_cur <-  c(-0.5,-0.002)

loglike <- loglike_sample(obs = 0, upper = upper_cur, 
                          logmu = logmu_cur, phi = phi_cur, gamma = gamma_cur)

loglike_upper1 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 1)


loglike_upper2 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 3)

loglike_upper3 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 5)
# Draw the traceplot
plot(loglike,type = "l", ylim = c(-1.5,0),
     xlab = "x_max", ylab = "loglikelihood",
     main = paste0("logmu = ",logmu_cur,
                   ", gamma_0 = ",gamma_cur[1],
                   ", gamma_1 = ",gamma_cur[2]))
abline(h = loglike_upper1, col = 2, lty = 2)
abline(h = loglike_upper2, col = 3, lty = 2)
abline(h = loglike_upper3, col = 4, lty = 2)
legend("topleft",legend = paste0("t = ",c(1,3,5)), lty =2, col = 2:4)

####################################
# Case 6: large mu and small gamma #
####################################
logmu_cur <- 5
upper_cur <- 1000
phi_cur <- 1
gamma_cur <-  c(-0.5,-0.002)

loglike <- loglike_sample(obs = 0, upper = upper_cur, 
                          logmu = logmu_cur, phi = phi_cur, gamma = gamma_cur)

loglike_upper1 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 1)


loglike_upper2 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 3)

loglike_upper3 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 5)
# Draw the traceplot
plot(loglike,type = "l", ylim = c(-5,0),
     xlab = "x_max", ylab = "loglikelihood",
     main = paste0("logmu = ",logmu_cur,
                   ", gamma_0 = ",gamma_cur[1],
                   ", gamma_1 = ",gamma_cur[2]))
abline(h = loglike_upper1, col = 2, lty = 2)
abline(h = loglike_upper2, col = 3, lty = 2)
abline(h = loglike_upper3, col = 4, lty = 2)
legend("topleft",legend = paste0("t = ",c(1,3,5)), lty =2, col = 2:4)

###################################
# Case 7: huge mu and large gamma #
###################################
logmu_cur <- 12
upper_cur <- 5000
phi_cur <- 1
gamma_cur <-  c(-0.5,-0.002)

loglike <- loglike_sample(obs = 0, upper = upper_cur, 
                          logmu = logmu_cur, phi = phi_cur, gamma = gamma_cur)

loglike_upper1 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 1)


loglike_upper2 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 3)

loglike_upper3 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 5)
# Draw the traceplot
plot(loglike,type = "l", ylim = c(-12,-6),
     xlab = "x_max", ylab = "loglikelihood",
     main = paste0("logmu = ",logmu_cur,
                   ", gamma_0 = ",gamma_cur[1],
                   ", gamma_1 = ",gamma_cur[2]))
abline(h = loglike_upper1, col = 2, lty = 2)
abline(h = loglike_upper2, col = 3, lty = 2)
abline(h = loglike_upper3, col = 4, lty = 2)
legend("topleft",legend = paste0("t = ",c(1,3,5)), lty =2, col = 2:4)


##########################################
# Case 8: moderate mu and moderate gamma #
##########################################
logmu_cur <- 1
upper_cur <- 100
phi_cur <- 1
gamma_cur <-  c(-0.5,-0.1)

loglike <- loglike_sample(obs = 0, upper = upper_cur, 
                          logmu = logmu_cur, phi = phi_cur, gamma = gamma_cur)

loglike_upper1 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 1)


loglike_upper2 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 3)

loglike_upper3 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 5)
# Draw the traceplot
plot(loglike,type = "l", ylim = c(-1.5,0),
     xlab = "x_max", ylab = "loglikelihood",
     main = paste0("logmu = ",logmu_cur,
                   ", gamma_0 = ",gamma_cur[1],
                   ", gamma_1 = ",gamma_cur[2]))
abline(h = loglike_upper1, col = 2, lty = 2)
abline(h = loglike_upper2, col = 3, lty = 2)
abline(h = loglike_upper3, col = 4, lty = 2)
legend("topleft",legend = paste0("t = ",c(1,3,5)), lty =2, col = 2:4)
abline(v = exp(logmu_cur),col = 2, lty = 3)

##########################################
# Case 8: moderate mu and moderate gamma #
##########################################
logmu_cur <- 10
upper_cur <- 100
phi_cur <- 1
gamma_cur <-  c(-0.5,-0.2)

loglike <- loglike_sample(obs = 0, upper = upper_cur, 
                          logmu = logmu_cur, phi = phi_cur, gamma = gamma_cur)

loglike_upper1 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 1)


loglike_upper2 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 3)

loglike_upper3 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 5)
# Draw the traceplot
plot(loglike,type = "l", ylim = c(-10,-5),
     xlab = "x_max", ylab = "loglikelihood",
     main = paste0("logmu = ",logmu_cur,
                   ", gamma_0 = ",gamma_cur[1],
                   ", gamma_1 = ",gamma_cur[2]))
abline(h = loglike_upper1, col = 2, lty = 2)
abline(h = loglike_upper2, col = 3, lty = 2)
abline(h = loglike_upper3, col = 4, lty = 2)
legend("topleft",legend = paste0("t = ",c(1,3,5)), lty =2, col = 2:4)
abline(v = 10,col = 2, lty = 3)

##########################################
# Case 9: moderate mu and moderate gamma #
##########################################
logmu_cur <- 0
upper_cur <- 100
phi_cur <- 1
gamma_cur <-  c(-0.5,-0.05)

loglike <- loglike_sample(obs = 0, upper = upper_cur, 
                          logmu = logmu_cur, phi = phi_cur, gamma = gamma_cur)

loglike_lower1 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 2)


loglike_upper2 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 3)


loglike_upper3 <- loglike_upper(obs = 0, logmu = logmu_cur, 
                                phi = phi_cur, gamma = gamma_cur, 
                                times = 5)
# Draw the traceplot
range(loglike)

plot(loglike,type = "l", ylim = c(-1,0),
     xlab = "x_max", ylab = "loglikelihood",
     main = paste0("logmu = ",logmu_cur,
                   ", gamma_0 = ",gamma_cur[1],
                   ", gamma_1 = ",gamma_cur[2]))

abline(h = (loglike_lower1 + loglike_upper2)/2, col = 2, lty = 2)
abline(h = loglike_upper2, col = 3, lty = 2)
abline(h = loglike_upper3, col = 4, lty = 2)
legend("topleft",legend = paste0("t = ",c(2.5,3,5)), lty =2, col = 2:4)
abline(v = 5,col = 2, lty = 3)
exp(logmu_cur + gamma_cur[2])

save.image("0509Upper_bound_loglike.RData")


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
