# Implement the benchmarked methods
library(ggplot2)
library(pROC)

rm(list=ls())
proj <- "simulation"
ver <- 1

########################
# Clustering accurarcy #
########################
methods <- c("DIFseq", "BUSseq", "Seurat", "Milo", "LIGER", "ZINBWaVE")
N_methods <- length(methods)
ARI_com <- rep(NA, length(methods))

for(i in 1:N_methods){
  
  method <- methods[i]
  ARI_com[i] <- unlist(read.table(paste0("Output/ARI_",method,"_simulation_v",ver,".txt")))
  
}


df_ARI <- data.frame(Method = factor(methods, levels = methods),
                     ARI = ARI_com)

pdf(file = paste0("Images/ARI_comparison_",proj,"_v",ver,"_adjusted.pdf"), 
    width = 10, height = 8)
p <- ggplot(data = df_ARI, mapping = aes(x = Method, fill = Method,
                                         y = (ARI - 0.4)/3 * 5)) +
  geom_bar(#color = "#939393", 
    stat = "identity"
    #position=position_dodge(
  ) +
  # scale_fill_manual(values=color_by_condition[2:3]) +
  scale_x_discrete(limits=rev) +
  scale_y_continuous(breaks = seq(0,1,1/6), labels = seq(0.4,1,0.1)) +
  coord_flip() +
  labs(x = "Method", y = "ARI"#,fill = "Method"
  ) +
  theme(axis.text=element_text(size=28),
        axis.title=element_text(size=32,face="bold"),
        legend.position="none")
p
dev.off()

#############
# ROC curve #
#############
# Load ROC results by MiloR
DE_ROC_DIFseq <- read.table(paste0("Output/DIFseq_cell_type_specific_DE_genes_for_ROC_v",ver,".txt"))

# Load ROC results by MiloR
DE_ROC_Milo <- read.table(paste0("Output/Milo_cell_type_specific_DE_genes_for_ROC_v",ver,".txt"))

# Load ROC results by Seurat
DE_ROC_Seurat <- read.table(paste0("Output/Seurat_cell_type_specific_DE_genes_for_ROC_v",ver,".txt"))

pdf(paste0("Images/ROC_DE_",proj,"_comparison_v",ver,".pdf"),width = 5, height = 5)
# ROC plot 1
roc(DE_ROC_DIFseq$DE, DE_ROC_DIFseq$logPrJ, plot=TRUE, legacy.axes=FALSE, percent=TRUE, col="#F8766D", lwd=2, print.auc=FALSE, cex.lab = 1.5, cex.axis = 1.5)

plot.roc(DE_ROC_Seurat$DE, DE_ROC_Seurat$pval, percent=TRUE, col="#00BA38", lwd=4, 
         print.auc=FALSE, add=TRUE, smooth = TRUE, lty = 2)

plot.roc(DE_ROC_Milo$DE, DE_ROC_Milo$p.adj, percent=TRUE, col="#619CFF", lwd=4, lty = 3, 
         print.auc=FALSE, add=TRUE, smooth = TRUE)
# title(main = "ROC comparison", line = 2.5)

# Add legend
legend("bottomright",
       legend=c("DIFseq (AUC: 100.0%)", "Seurat (AUC: 92.6%)", "Milo (AUC: 88.6%)"),
       lty = c(1,2,3),
       col=c("#F8766D", "#00BA38", "#619CFF"),
       lwd=2, xpd = TRUE, bty = "n")
dev.off()