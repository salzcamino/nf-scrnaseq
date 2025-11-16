process CLUSTERING {
    tag "clustering"
    label 'process_medium'
    publishDir "${params.outdir}/clustering", mode: params.publish_dir_mode

    input:
    path adata
    val run_leiden
    val run_louvain
    val leiden_resolution
    val louvain_resolution
    val cluster_key

    output:
    path "clustered.h5ad", emit: adata
    path "cluster_assignments.csv", emit: clusters
    path "clustering_plots.pdf", emit: plots
    path "clustering_summary.txt", emit: summary

    shell:
    '''
#!/usr/bin/env python3
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

print("Loading data with dimensionality reductions...")
adata = sc.read_h5ad('!{adata}')

n_cells = adata.n_obs
n_genes = adata.n_vars
print(f"Input: {n_cells} cells, {n_genes} genes")

# Parse parameters
run_leiden_bool = '!{run_leiden}'.lower() == 'true'
run_louvain_bool = '!{run_louvain}'.lower() == 'true'
leiden_res = float('!{leiden_resolution}')
louvain_res = float('!{louvain_resolution}')
cluster_key = '!{cluster_key}'

results = {}

# Check that neighbor graph exists
if 'neighbors' not in adata.uns:
    print("ERROR: Neighbor graph not found. Run dimensionality reduction first.")
    raise ValueError("Neighbor graph not computed")

# Leiden clustering
if run_leiden_bool:
    print(f"Running Leiden clustering (resolution={leiden_res})...")
    sc.tl.leiden(adata, resolution=leiden_res, key_added='leiden')
    n_leiden_clusters = len(adata.obs['leiden'].unique())
    results['leiden'] = {
        'n_clusters': n_leiden_clusters,
        'resolution': leiden_res
    }
    print(f"Leiden: Found {n_leiden_clusters} clusters")

# Louvain clustering
if run_louvain_bool:
    print(f"Running Louvain clustering (resolution={louvain_res})...")
    sc.tl.louvain(adata, resolution=louvain_res, key_added='louvain')
    n_louvain_clusters = len(adata.obs['louvain'].unique())
    results['louvain'] = {
        'n_clusters': n_louvain_clusters,
        'resolution': louvain_res
    }
    print(f"Louvain: Found {n_louvain_clusters} clusters")

# Set default cluster key if not specified
if cluster_key == 'auto':
    if run_leiden_bool:
        cluster_key = 'leiden'
    elif run_louvain_bool:
        cluster_key = 'louvain'
    else:
        cluster_key = None

# Save cluster assignments
print("Saving cluster assignments...")
cluster_df = pd.DataFrame(index=adata.obs_names)
if run_leiden_bool:
    cluster_df['leiden'] = adata.obs['leiden'].values
if run_louvain_bool:
    cluster_df['louvain'] = adata.obs['louvain'].values
cluster_df.to_csv('cluster_assignments.csv')

# Create plots
print("Generating clustering plots...")

# Determine number of plots needed
n_cluster_methods = int(run_leiden_bool) + int(run_louvain_bool)
has_umap = 'X_umap' in adata.obsm
has_pca = 'X_pca' in adata.obsm

# Create figure with appropriate size
n_rows = max(2, n_cluster_methods)
n_cols = 3 if has_umap else 2
fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
if n_rows == 1:
    axes = axes.reshape(1, -1)

plot_row = 0

# Leiden plots
if run_leiden_bool:
    # Cluster distribution
    ax = axes[plot_row, 0]
    cluster_counts = adata.obs['leiden'].value_counts().sort_index()
    ax.bar(range(len(cluster_counts)), cluster_counts.values, color='steelblue', alpha=0.7)
    ax.set_xlabel('Leiden Cluster')
    ax.set_ylabel('Number of cells')
    ax.set_title(f'Leiden Clustering (n={n_leiden_clusters})')
    ax.set_xticks(range(len(cluster_counts)))
    ax.set_xticklabels(cluster_counts.index, rotation=45)

    # PCA with clusters
    if has_pca:
        ax = axes[plot_row, 1]
        sc.pl.pca(adata, color='leiden', ax=ax, show=False, title='Leiden clusters (PCA)', legend_loc='right margin')

    # UMAP with clusters
    if has_umap:
        ax = axes[plot_row, 2] if has_pca else axes[plot_row, 1]
        sc.pl.umap(adata, color='leiden', ax=ax, show=False, title='Leiden clusters (UMAP)', legend_loc='right margin')

    plot_row += 1

# Louvain plots
if run_louvain_bool:
    # Cluster distribution
    ax = axes[plot_row, 0]
    cluster_counts = adata.obs['louvain'].value_counts().sort_index()
    ax.bar(range(len(cluster_counts)), cluster_counts.values, color='darkorange', alpha=0.7)
    ax.set_xlabel('Louvain Cluster')
    ax.set_ylabel('Number of cells')
    ax.set_title(f'Louvain Clustering (n={n_louvain_clusters})')
    ax.set_xticks(range(len(cluster_counts)))
    ax.set_xticklabels(cluster_counts.index, rotation=45)

    # PCA with clusters
    if has_pca:
        ax = axes[plot_row, 1]
        sc.pl.pca(adata, color='louvain', ax=ax, show=False, title='Louvain clusters (PCA)', legend_loc='right margin')

    # UMAP with clusters
    if has_umap:
        ax = axes[plot_row, 2] if has_pca else axes[plot_row, 1]
        sc.pl.umap(adata, color='louvain', ax=ax, show=False, title='Louvain clusters (UMAP)', legend_loc='right margin')

    plot_row += 1

# Hide unused axes
for i in range(plot_row, n_rows):
    for j in range(n_cols):
        axes[i, j].axis('off')

plt.suptitle('Clustering Results', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('clustering_plots.pdf', dpi=150, bbox_inches='tight')
plt.close()

# Save clustered data
print("Saving clustered data...")
adata.write('clustered.h5ad')

# Write summary
summary = f"""Clustering Summary
==================
Input cells: {n_cells}
Input genes: {n_genes}

Parameters:
  Run Leiden: {run_leiden_bool}
  Run Louvain: {run_louvain_bool}
  Leiden resolution: {leiden_res}
  Louvain resolution: {louvain_res}
  Primary cluster key: {cluster_key}

"""

if run_leiden_bool:
    summary += f"""Leiden Clustering Results:
  Number of clusters: {results['leiden']['n_clusters']}
  Resolution: {results['leiden']['resolution']}

  Cluster sizes:
"""
    for cluster_id in sorted(adata.obs['leiden'].unique(), key=lambda x: int(x)):
        count = (adata.obs['leiden'] == cluster_id).sum()
        pct = 100 * count / n_cells
        summary += f"    Cluster {cluster_id}: {count} cells ({pct:.1f}%)\\n"
    summary += "\\n"

if run_louvain_bool:
    summary += f"""Louvain Clustering Results:
  Number of clusters: {results['louvain']['n_clusters']}
  Resolution: {results['louvain']['resolution']}

  Cluster sizes:
"""
    for cluster_id in sorted(adata.obs['louvain'].unique(), key=lambda x: int(x)):
        count = (adata.obs['louvain'] == cluster_id).sum()
        pct = 100 * count / n_cells
        summary += f"    Cluster {cluster_id}: {count} cells ({pct:.1f}%)\\n"
    summary += "\\n"

summary += f"""
Data saved to: clustered.h5ad
Cluster assignments saved to: cluster_assignments.csv

Available cluster annotations:
"""
if run_leiden_bool:
    summary += "  - adata.obs['leiden']\\n"
if run_louvain_bool:
    summary += "  - adata.obs['louvain']\\n"

with open('clustering_summary.txt', 'w') as f:
    f.write(summary)

print(summary)
print("Clustering complete!")
    '''
}
