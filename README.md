# Differential inference for scRNA-seq data implementation

## Overview

As comparisons between biological conditions or treatments are one of the most important methodologies for scientific discoveries, single-cell RNA-sequencing (scRNA-seq) experiments with cells collected from multiple conditions have become prevalent nowadays. However, rigorous statistical methods for comparing scRNA-seq data collected from multiple conditions are lacking. 

Therefore, we propose an interpretable hierarchical Bayesian model **DIFseq** that is able to simultaneously correct batch effects, perform cell type clustering, impute dropout events, and estimate batch effects, cell type effects and cell-type-specific condition effects. Via joint modeling, DIFseq accounts for all the uncertainties in the analysis and enables rigorous hypothesis testing. By borrowing ideas from deep learning, we implement a Monte Carlo Expectation-Maximization (MCEM) algorithm with stochastic updates so that the algorithm is scalable to hundreds of thousands of cells.

This repository introduces how to reproduce the results in our manuscript. First, we should compile the source code of DIFseq and then reproduce the figures and tables for three simulation settings and two real datasets as the following instructions.  

## Compiling the source code

The algorithm relies on GNU scietific library. To compile the source code on your computer, please modify the directory of GSL in the *Makefile* file. Then, you can compile the C++ source code by

```
make
```

We apply [OMPRNG](https://homepage.divms.uiowa.edu/~mbognar/omprng/) to generate random number in parallel.


## Simulation

Our manuscript includes three simulation settings, including

   - Setting 1: Assume that no batch effect exists between different batches and the cell-size factors of all cells are the same as if the batch effects correction and cell-wise normalization were perfect. 
   - Setting 2: Incorperate batch effects and allow different cell-size factors
   - Setting 3: In the first two settings, all conditions exist in all batches. Setting 3 allows a flexible experimental design where some conditions are missing in some batches.

Folder `vX` contains the reproduction code of setting `X`, where `X` takes the values of 1, 2 and 3 corresponding to three different settings. 

### Generate synthetic datasets

To generate the synthetic datasets of three settings, please run `Simulate_data_simulation_vX.R` to generate the synthetic datasets together with dimension information and metadata in the `RawCountData` folder.

### Analysis by DIFseq 

Then, run `run_DIFseq_simulation.sh` to implement the MCEM algorithm and conduct statistical inference. The parameter estimation will be stored in the `Inferece_KY` folder, where `Y` denotes the number of cell types. Here, we vary the number of cell types from 3 to 8 for each settings. 

*Note*: Each line in the `run_DIFseq_simulation.sh` corresponds to one cell type number, so you can also run different lines separately or only run one single line (like the optimal cell type number).

After finishing the MCEM algorithm, `DIFseq_analysis.R` will

   - Draw the line chart of BICs to select the optimal number of cell types;
   - Draw the bar plot of the true and estimated cell-type proprotions;
   - Draw the heatmaps of the true and estimated cell type effects, batch effects and condition effects;
   - Draw the heatmaps of raw read counts and corrected read counts;
   - Draw the UMAP of corrected read counts colored by cell type, by batch and by condition, respecitively;
   - Compute the test statistics and p-values for differential abundance analysis;
   - Draw the UMAP of corrected read counts colored by negative log p-values to visualize the differential abundance results;
   - Conduct differnetial expression analysis to generate the ROC curve of cell-type-specific DE genes.

### Benchmarked methods

To evaluate the **clustering performance**, we benmark DIFseq with some state-of-the-arts methods. The corresponding R code of each method is listed as below:

   - BUSseq: `BUSseq_clustering_analysis.R`;
   - Seurat: `Seurat_DE.R` for v1 and `Batch_effects_correction.R` for v2 and v3; 
   - Milo: `Milo_DA_analysis.R`;
   - liger: `LIGER_clustering_analysis.R`;
   - ZINB-WaVE: `ZINBWaVE_clustering_analysis.R`.
   
 
We also benchmark the **differential analysis (DA) analysis** of DIFseq with the state-of-the-arts DA methods. The implementation of each method is listed as below:

   - DA-seq is only applicable for the first two simulation settings by running `DAseq_DA_analysis.R`;
   - MELD is also applicable for the frist two settings by first running `MELD_DA_analysis.py` and then running `MELD_DA_analysis.R` to draw the UMAP of relative likelihood;
   - Milo works for all the three settings by running `Milo_DA_analysis.R`.

*Note:* Before conducting DA analysis on settings 2 and 3, you need to run `Batch_effects_correction.R` to remove batch effects and generate the corrected PCs by Seurat or MNN. In comparison, settting 1 is free of batch effects. You can conduct DA analysis directly. 

To evaluate the **differential expression (DE)** analysis，we compare the ROC of identified cell-type-specific DE genes by DIFseq, Seurat and Milo, respectively. We also compare the area under the ROC curves (AUC) across three methods. The code of these three methods also involves the **differential expression** steps.

Finally, `Method_comparison.R` covers

   - Drawing the bar plot of ARI in cell type clustering;
   - Plotting the combined ROC curves of DIFseq, Seurat and Milo in differential expression analysis.

## Real data analysis

### Download raw count data 

Please download the raw count data of the human pancreas study and the human COVID-19 study from [the OneDrive link](https://cuhko365-my.sharepoint.com/:f:/g/personal/songfangda_cuhk_edu_cn/ErmSyLVRFhlDvNxETCSCi8gB5So6RqxJ9KH4R7B99s0UiQ?e=eWaCRX) with password `DIFseq2024`. To run the reproduction code, please place the `RawCountData` folders into their corresponding `Pancreas` and `Covid` folders, respectively.

### Pancreas study

After placing the `RawCountData` folder inside the `Pancreas` folder, please run `run_DIFseq_pancreas.sh` to implement the MCEM algorithm and conduct statistical inference for different cell type numbers. 

To analyze the pancreas data by DIFseq, you can run `DIFseq_analysis.R` to

   - Draw the scatter plot of BICs for different cell type numbers to select the optimal number of cell types;
   - Draw the bar plot of estimated cell type proportions of each donor;
   - Conduct differential abundance analysis for the two batches with two conditions;
   - Conduct differential expression analysis and draw the condition effects of cell-type-specific DE genes;
   - Generate the corrected read count data and draw dot plot for the corrected expression levels of marker genes;
   - Draw violion plots of the identified DE genes for Beta cells.

Then, we also conduct differential expression analysis by DA-seq (`DAseq_DA_analysis.R`), MELD (`MELD_DA_analysis.py` and `MELD_DA_analysis.R`) and Milo (`Milo_DA_analysis.R`) to the last two batches with two conditions.

### COVID-19 study

After placing the `RawCountData` folder into the `Covid` folder, please run `run_DIFseq_covid.sh` to implement the MCEM algorithm and conduct statistical inference for different cell type numbers. 

To analyze the covid-19 data by DIFseq, you can run `DIFseq_analysis.R` to

   - Draw the scatter plot of BICs for different cell type numbers to select the optimal number of cell types;
   - Annotate the identified cell types by their average anti-body (ADT) levels and draw the heatmap of the scaled ADT levels;
   - Draw the bar plot of the estimated cell type proportions;
   - Conduct differential expression analysis and draw the condition effects of cell-type-specific DE genes as well as the bar plot of DE gene numbers for two conditions;
   - Generate the corrected read count data and draw the dot plot for the corrected expression levels of immune-related genes;
   - Conduct differential abundance analysis between any pair of three conditions.

Then, we also apply Milo (`Milo_DA_analysis.R`) and MELD (`MELD_DA_analysis.R`) to conduct differential abundance analysis.





