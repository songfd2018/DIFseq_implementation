import pandas as pd
import numpy as np
import graphtools as gt
import matplotlib.pyplot as plt
import cmocean
import magic
import phate
import scprep
import meld
import sklearn
import tempfile
import os
import scanpy as sc
import umap

# Make sure plots & clusters reproducible
np.random.seed(42)

########################
# Load PC and metadata #
########################
print('Loading count data and metadata...')
proj = 'covid'
ver = 1

# Load from csv
## PCA 
pca_file = 'RawCountData/PCA_MNN_{proj}_v{ver}.csv'
data_pca = scprep.io.load_csv(pca_file.format(proj = proj, ver = ver),
                               cell_names=True, gene_names=True, sparse=False)


## UMAP 
umap_file = 'RawCountData/UMAP_MNN_{proj}_v{ver}.csv'
df_umap = scprep.io.load_csv(umap_file.format(proj = proj, ver = ver),
                               cell_names=True, gene_names=True, sparse=False)
data_umap = df_umap.to_numpy()


## Metedata
meta_file = 'RawCountData/metadata_{proj}_v{ver}.csv'
metadata = scprep.io.load_csv(meta_file.format(proj = proj, ver = ver),
                               cell_names=True, gene_names=True, sparse=False)

############
# Run MELD #
############
# beta is the amount of smoothing to do for density estimation
# knn is the number of neighbors used to set the kernel bandwidth
meld_op = meld.MELD(beta=67, knn=7)
metadata['Pair_ID'] = ['_'.join(i) for i in metadata[['Severity','Batch']].values]

print('The first several row of metadata is given by:')
print(metadata.head())

sample_densities = meld_op.fit_transform(data_pca, sample_labels=metadata['Pair_ID'])

def replicate_normalize_densities(sample_densities, replicate):
    replicates = np.unique(replicate)
    sample_likelihoods = sample_densities.copy()
    for rep in replicates:
        curr_cols = sample_densities.columns[[col.endswith(rep) for col in sample_densities.columns]]
        print(curr_cols)
        sample_likelihoods[curr_cols] = sklearn.preprocessing.normalize(sample_densities[curr_cols], norm='l1')
    return sample_likelihoods


# print(metadata.head())
print('The first several row of sample_densities is given by:')
print(sample_densities.head())
# print(np.unique(metadata['batch']))


sample_likelihoods = replicate_normalize_densities(sample_densities, metadata['Batch'])
# print(type(sample_likelihoods))

print('The first several row of sample_likelihoods is given by:')
print(sample_likelihoods.head())

print('The dimension of sample_likelihoods is given by:')
print(sample_likelihoods.shape)

# Store the relative sample likelihood to a csv file
os.makedirs('Output', exist_ok=True)

sample_like_file = 'Output/sample_like_MELD_{proj}_v{ver}.csv'
sample_likelihoods.to_csv(sample_like_file.format(proj = proj, ver = ver))  

##################################################
# Draw the relative likelihood within each batch #
##################################################
umap_plot = 'Images/UMAP_colored_by_{group}_{proj}_v{ver}.jpg'
experimental_samples = ['Healthy_Batch1', 'Healthy_Batch2']

metadata['Healthy_likelihood'] = sample_likelihoods[experimental_samples].mean(axis=1).values

experimental_samples = ['Moderate_Batch1', 'Moderate_Batch2']

metadata['Moderate_likelihood'] = sample_likelihoods[experimental_samples].mean(axis=1).values

experimental_samples = ['Severe_Batch1', 'Severe_Batch2']

metadata['Severe_likelihood'] = sample_likelihoods[experimental_samples].mean(axis=1).values

fig, axes = plt.subplots(1,2, figsize=(13,6))
for i, ax in enumerate(axes):
    curr_sample = experimental_samples[i]

    # print(i)
    # print(curr_sample)
    # print(data_umap[:5,:])
    # print(meld.get_meld_cmap())

    scprep.plot.scatter2d(data_umap, c=sample_likelihoods[curr_sample], cmap=meld.get_meld_cmap(),
                          vmin=0, vmax=1,
                          title=curr_sample, ticks=False, ax=ax)

fig.savefig(umap_plot.format(group = 'likelihood_healthy', proj = proj, ver = ver))

## Examining the distribution of chd likelihood values in published clusters
# metadata['Cluster_l1'] = scprep.utils.sort_clusters_by_values(metadata['Cluster_l1'], metadata['Chd_likelihood'])

meta_out_file = 'Output/metadata_MELD_{proj}_v{ver}.csv'
metadata.to_csv(meta_out_file.format(proj = proj, ver = ver))  

# Draw jitter plot for all the cells
sample_cmap = {'Healthy_Batch1': '#3182bd',
               'Healthy_Batch2': '#08519c',
               'Moderate_Batch1': '#fec44f',
               'Moderate_Batch2': '#feb24c',
               'Severe_Batch1': '#fb6a4a',
               'Severe_Batch2': '#de2d26'}

jitter_plot = 'Images/MELD_jitter_plot_{proj}_v{ver}.jpg'


print('The first several row of metadata is given by:')
print(metadata.head())

# metadata['Cond_index'] = [1 if sl.startswith('N') else 0 for sl in metadata['Condition']]

fig, axes = plt.subplots(1, 3, figsize=(30,10))

metadata['Cluster_l1'] = scprep.utils.sort_clusters_by_values(metadata['Cluster_l1'], metadata['Healthy_likelihood'])
# See example usage: https://scprep.readthedocs.io/en/stable/examples/jitter.html
scprep.plot.jitter(metadata['Cluster_l1'], metadata['Healthy_likelihood'], c=metadata['Pair_ID'], 
                   cmap=sample_cmap,legend=False, plot_means=False, xlabel=False, ylabel='Mean healthy likelihood',
                   ax=axes[0])

### This code will plot the ratio of tyr:chd cells per cluster
metadata['Cond_index'] = [1 if sl.startswith('H') else 0 for sl in metadata['Severity']]
means = metadata.groupby('Cluster_l1')['Cond_index'].mean()
axes[0].scatter(means.index, means - np.mean(metadata['Cond_index']) + 0.5, color='#e0d8f0', edgecolor='k', s=100)
# axes[0].scatter(means.index, means, color='#e0d8f0', edgecolor='k', s=100)

print("The first subgraph")
print(metadata.head())
print(means)

# Axis tick labels
axes[0].set_xticklabels(metadata.set_index('Cluster_l1')['Seurat_Level1'].drop_duplicates().sort_index(), rotation=90)
axes[0].set_ylim(0,1)

metadata['Cluster_l1'] = scprep.utils.sort_clusters_by_values(metadata['Cluster_l1'], metadata['Moderate_likelihood'])
scprep.plot.jitter(metadata['Cluster_l1'], metadata['Moderate_likelihood'], c=metadata['Pair_ID'], 
                   cmap=sample_cmap,legend=False, plot_means=False, xlabel=False, ylabel='Mean moderate likelihood',
                   ax=axes[1])

### This code will plot the ratio of tyr:chd cells per cluster
metadata['Cond_index'] = [1 if sl.startswith('M') else 0 for sl in metadata['Severity']]
means = metadata.groupby('Cluster_l1')['Cond_index'].mean()
axes[1].scatter(means.index, means - np.mean(metadata['Cond_index']) + 0.5, color='#e0d8f0', edgecolor='k', s=100)
# axes[1].scatter(means.index, means, color='#e0d8f0', edgecolor='k', s=100)


# Axis tick labels
axes[1].set_xticklabels(metadata.set_index('Cluster_l1')['Seurat_Level1'].drop_duplicates().sort_index(), rotation=90)
axes[1].set_ylim(0,1)

metadata['Cluster_l1'] = scprep.utils.sort_clusters_by_values(metadata['Cluster_l1'], metadata['Severe_likelihood'])
scprep.plot.jitter(metadata['Cluster_l1'], metadata['Severe_likelihood'], c=metadata['Pair_ID'], 
                   cmap=sample_cmap,legend=False, plot_means=False, xlabel=False, ylabel='Mean severe likelihood',
                   ax=axes[2])

### This code will plot the ratio of tyr:chd cells per cluster
metadata['Cond_index'] = [1 if sl.startswith('S') else 0 for sl in metadata['Severity']]
means = metadata.groupby('Cluster_l1')['Cond_index'].mean()
axes[2].scatter(means.index, means - np.mean(metadata['Cond_index']) + 0.5, color='#e0d8f0', edgecolor='k', s=100)
# axes[2].scatter(means.index, means, color='#e0d8f0', edgecolor='k', s=100)

print("The third subgraph")
print(means)
print(np.mean(metadata['Cond_index']))

# Axis tick labels
axes[2].set_xticklabels(metadata.set_index('Cluster_l1')['Seurat_Level1'].drop_duplicates().sort_index(), rotation=90)
axes[2].set_ylim(0,1)


fig.savefig(jitter_plot.format(proj = proj, ver = ver))
