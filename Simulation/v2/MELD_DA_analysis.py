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
proj = 'simulation'
ver = 2

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

# print(type(counts_file))
# print(type(data))

## Metedata
meta_file = 'RawCountData/metadata_{proj}_v{ver}.csv'
metadata = scprep.io.load_csv(meta_file.format(proj = proj, ver = ver),
                               cell_names=True, gene_names=True, sparse=False)

# print(type(metadata))
# print(metadata.head())

############
# Run MELD #
############
# beta is the amount of smoothing to do for density estimation
# knn is the number of neighbors used to set the kernel bandwidth
meld_op = meld.MELD(beta=67, knn=7)
metadata['Pair_ID'] = ['_'.join(i) for i in metadata[['Condition','Batch']].values]

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
fig, axes = plt.subplots(1,4, figsize=(13,4))

umap_plot = 'Images/UMAP_colored_by_{group}_{proj}_v{ver}.jpg'
experimental_samples = ['Condition_1_Batch_1', 'Condition_1_Batch_2', 'Condition_1_Batch_3', 'Condition_1_Batch_4']

for i, ax in enumerate(axes):
    curr_sample = experimental_samples[i]

    # print(i)
    # print(curr_sample)
    # print(data_umap[:5,:])
    # print(meld.get_meld_cmap())

    scprep.plot.scatter2d(data_umap, c=sample_likelihoods[curr_sample], cmap=meld.get_meld_cmap(),
                          vmin=0, vmax=1,
                          title=curr_sample, ticks=False, ax=ax)

fig.savefig(umap_plot.format(group = 'likelihood', proj = proj, ver = ver))

# They use the average likelihood of the chordin samples as the measure of the perturbation
fig, axes = plt.subplots(1,2, figsize=(8.7,4))

scprep.plot.scatter2d(data_umap, c=sample_likelihoods[experimental_samples].mean(axis=1), 
                      cmap=meld.get_meld_cmap(), vmin=0, vmax=1,
                      title='Mean', ticks=False, ax=axes[0])
scprep.plot.scatter2d(data_umap, c=sample_likelihoods[experimental_samples].std(axis=1), vmin=0, 
                      cmap='inferno', title='St. Dev.', ticks=False, ax=axes[1])

fig.tight_layout()

metadata['Chd_likelihood'] = sample_likelihoods[experimental_samples].mean(axis=1).values


## Examining the distribution of chd likelihood values in published clusters
metadata['CellType'] = scprep.utils.sort_clusters_by_values(metadata['CellType'], metadata['Chd_likelihood'])

meta_out_file = 'Output/metadata_MELD_{proj}_v{ver}.csv'
metadata.to_csv(meta_out_file.format(proj = proj, ver = ver))  