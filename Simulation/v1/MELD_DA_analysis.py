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
ver = 1

# Load from csv
## PCA 
pca_file = 'RawCountData/PCA_corrected_{proj}_v{ver}.csv'
data_pca = scprep.io.load_csv(pca_file.format(proj = proj, ver = ver),
                               cell_names=True, gene_names=True, sparse=False)


## UMAP 
umap_file = 'RawCountData/UMAP_corrected_{proj}_v{ver}.csv'
df_umap = scprep.io.load_csv(umap_file.format(proj = proj, ver = ver),
                               cell_names=True, gene_names=True, sparse=False)
data_umap = df_umap.to_numpy()

# print(type(counts_file))
# print(type(data))

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
metadata['Pair_ID'] = ['_'.join(i) for i in metadata[['Condition','Batch']].values]

print('The first several row of metadata is given by:')
print(metadata.head())

sample_densities = meld_op.fit_transform(data_pca, sample_labels=metadata['Pair_ID'])

# print(type(sample_densities))
# print(sample_densities.shape)
# print(sample_densities.head())
# print(sample_densities.columns)

# sample_densities 
# Row for cells
# Column for samples
# Normalize the densities of samples belonging to the same group
# How to understand the replicate? Pairwise studies? Confirm from the paper

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

meta_out_file = 'Output/metadata_MELD_{proj}_v{ver}.csv'
metadata.to_csv(meta_out_file.format(proj = proj, ver = ver))  


## Sort cells by cluster and Chd_likelihood
# metadata['CellType'] = scprep.utils.sort_clusters_by_values(metadata['CellType'], metadata['Chd_likelihood'])

fig, ax = plt.subplots(1, figsize=(10,10))

# See example usage: https://scprep.readthedocs.io/en/stable/examples/jitter.html
scprep.plot.jitter(metadata['CellType'], metadata['Chd_likelihood'], c=metadata['Sample'], 
                   cmap=sample_cmap,legend=False, plot_means=False, xlabel=False, ylabel='Mean chd likelihood',
                   ax=ax)

### This code will plot the ratio of tyr:chd cells per cluster
means = metadata.groupby('CellType')['Condition'].mean()
ax.scatter(means.index, means - np.mean(metadata['Condition']) + 0.5, color='#e0d8f0', edgecolor='k', s=100)

fig.tight_layout()

###############
# Perform VFC #
###############
data_cluster_UMAP = {}

for cluster in clusters:
    data_cluster_UMAP[cluster] = data_umap.loc[metadata['CellType'] == cluster, :] 


fig,axes= plt.subplots(5,4, figsize=(4*3, 5*3))

for i , ax in enumerate(axes.flatten()):
    if not i < len(clusters):
        ax.axis('off')
        continue
    curr_cluster = clusters[i]
    curr_UMAP = data_cluster_UMAP[curr_cluster]
    
    scprep.plot.scatter2d(curr_UMAP, 
                          c=metadata['Chd_likelihood'].loc[metadata['CellType'] == curr_cluster], 
                          cmap=meld.get_meld_cmap(), vmin=0, vmax=1,
                         ax=ax, ticks=False, 
                          title='Cluster {} ({})'.format(curr_cluster, curr_UMAP.shape[0]), 
                          legend=False, fontsize=10)

fig.savefig(umap_plot.format(group = 'Cluster', proj = proj, ver = ver))

# Build a VFC operator for each cluster
np.random.seed(0)
vfc_op_per_cluster = {}

for cluster in np.unique(metadata['CellType']):
    curr_G = gt.Graph(data_pca.loc[metadata['CellType'] == cluster], use_pygsp=True)
    curr_G.compute_fourier_basis()
    curr_sample_labels = metadata['Condition'].loc[metadata['CellType'] == cluster]
    curr_likelihood = metadata['Chd_likelihood'].loc[metadata['CellType'] == cluster]
    curr_vfc = meld.VertexFrequencyCluster(n_clusters = 3)
    curr_vfc.fit_transform(curr_G, curr_sample_labels, curr_likelihood)
    vfc_op_per_cluster[cluster] = curr_vfc

subclustering_results = {}
for cluster in np.unique(metadata['CellType']):
    curr_vfc = vfc_op_per_cluster[cluster]
    clusters_by_n = {}
    for n in [2,3,4,5]:
        clusters_by_n[n] = curr_vfc.predict(n)
    subclustering_results[cluster] = clusters_by_n


fig, axes= plt.subplots(19,5, figsize=(15, 19*3))

for i, r_ax in enumerate(axes):
    cluster = clusters[i]
    curr_UMAP = data_cluster_UMAP[cluster]

    for i, n in enumerate([0,2,3,4,5]):
        if i == 0:
            cvec = metadata['chd_likelihood'].loc[metadata['clusterID'] == cluster]
            ylabel= str(cluster) + ' ({})'.format(len(cvec))

            cmap= meld.get_meld_cmap()
            vmin=0
            vmax=1           
        else:
            cvec = subclustering_results[cluster][n]
            cmap = ylabel = vmin = vmax = None
            
        if np.sum(metadata['clusterID'] == cluster) < (metadata.shape[0] * 0.01):
            plt.setp(r_ax[i].spines.values(), color='lightgrey', linewidth=4)
        scprep.plot.scatter2d(curr_phate, c=cvec, cmap=cmap, vmin=vmin, vmax=vmax,
                         ax=r_ax[i], ticks=False, ylabel=ylabel, fontsize=8, legend=False)
        # Red outline for clusters that we decide to subcluster
        if cluster in picked_clusters:
            if picked_clusters[cluster] == n:
                plt.setp(r_ax[i].spines.values(), color='r', linewidth=4)

# Load corrected data
count_data_file = '../RawCountData/count_data_{proj}_v{ver}.csv'
count_data = scprep.io.load_csv(count_data_file.format(proj = proj, ver = ver),
                               cell_names=True, gene_names=True, sparse=False)

# celltype_marker_genes = {}
DE_out_file = 'Output/DE_MELD_{proj}_v{ver}_celltype{type}.csv'

# Iterate over VFC subclusters
for curr_celltype in np.unique(metadata['CellType']):

    celltype_data = count_data.loc[metadata['CellType'] == curr_celltype, :]
    # Create mask for conditions
    celltype_condition = metadata['Condition'].loc[metadata['CellType'] == curr_celltype]
    
    # Rank test comparing each subcluster to all other cells in cluster
    curr_de_results = de.test.two_sample(data = celltype_data, 
                                         # gene_names = data.columns, 
                                         grouping = celltype_condition,
                                         test='rank').summary()
                            
    curr_de_results.to_csv(DE_out_file.format(proj = proj, ver = ver, type = curr_celltype))
