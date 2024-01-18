######################################
# Apply Milo to the Covid-19 dataset #
######################################
library(dplyr)
library(patchwork)
library(scater)
library(miloR)

rm(list=ls())
set.seed(1234)
proj <- "covid"
ver <- 1

########################################
# load the raw count data and metadata #
########################################
load(paste0("Output/MNN_analysis_",proj,"_v",ver,".RData"))

##################################
# Differential abundance testing #
##################################
# Create a Milo object
covid_milo <- Milo(fastCorrect)
covid_milo

# Construct KNN graph
covid_milo <- buildGraph(covid_milo, k = 30, d = 30, reduced.dim = "pca.corrected")
covid_milo <- makeNhoods(covid_milo, prop = 0.1, k = 30, d=30, refined = TRUE, reduced_dims = "pca.corrected")

# Counting cells in neighbourhoods
covid_milo <- countCells(covid_milo, meta.data = metadata, samples="Pair")

covid_design <- metadata[,c("Pair", "Severity", "Batch")]
## Convert batch info from integer to factor
covid_design$Batch <- factor(covid_design$Batch)
covid_design$Severity <- factor(covid_design$Severity, levels = c("Healthy", "Moderate", "Severe"))
covid_design <- distinct(covid_design)
rownames(covid_design) <- covid_design$Pair

# 3.6 computing neighbourhood connectivity
covid_milo <- calcNhoodDistance(covid_milo, d=30, reduced.dim = "pca.corrected")

# 3.7 testing
da_results_allcom <- testNhoods(covid_milo, design = ~ Severity,
                                model.contrasts = c("SeverityModerate",
                                                    "SeveritySevere",
                                                    "SeverityModerate-SeveritySevere"),
                                design.df = covid_design, reduced.dim = "UMAP")

head(da_results_allcom)

colnames(da_results_allcom)[1:3] <- paste0("logFC_",c("Healthy_Moderate", "Healthy_Severe", "Moderate_Severe"))

# Output the number of identified DA neighbors
sum(da_results_allcom$SpatialFDR < 0.1)
