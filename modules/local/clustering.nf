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
    val run_seurat
    val run_celda
    val seurat_resolution
    val celda_L
    val celda_K

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

# R-based clustering parameters
run_seurat_bool = '!{run_seurat}'.lower() == 'true'
run_celda_bool = '!{run_celda}'.lower() == 'true'
seurat_res = float('!{seurat_resolution}')
celda_L_str = '!{celda_L}'
celda_K_str = '!{celda_K}'

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

# Seurat SNN clustering (R-based)
if run_seurat_bool:
    print(f"Running Seurat SNN clustering (resolution={seurat_res})...")
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.packages import importr

        pandas2ri.activate()

        # Load Seurat
        seurat = importr('Seurat')

        # Create Seurat object from count matrix
        # Get raw counts if available
        if 'counts' in adata.layers:
            counts = adata.layers['counts']
        else:
            counts = adata.X

        # Convert to dense if sparse
        if hasattr(counts, 'toarray'):
            counts_dense = counts.toarray()
        else:
            counts_dense = counts

        # Create R matrix (genes x cells)
        counts_df = pd.DataFrame(
            counts_dense.T,
            index=adata.var_names,
            columns=adata.obs_names
        )

        ro.r("""
        run_seurat_clustering <- function(counts_df, resolution, n_pcs) {
            # Create Seurat object
            seurat_obj <- CreateSeuratObject(counts = as.matrix(counts_df))

            # Normalize
            seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

            # Find variable features
            seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)

            # Scale
            seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)

            # PCA
            seurat_obj <- RunPCA(seurat_obj, npcs = n_pcs, verbose = FALSE)

            # Find neighbors
            seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_pcs, verbose = FALSE)

            # Find clusters
            seurat_obj <- FindClusters(seurat_obj, resolution = resolution, verbose = FALSE)

            # Return cluster assignments
            return(as.character(Idents(seurat_obj)))
        }
        """)

        n_pcs_for_seurat = min(30, n_cells - 1)
        seurat_clusters = ro.r['run_seurat_clustering'](
            pandas2ri.py2rpy(counts_df),
            seurat_res,
            n_pcs_for_seurat
        )

        # Convert to pandas and add to adata
        adata.obs['seurat_clusters'] = pd.Categorical(list(seurat_clusters))
        n_seurat_clusters = len(adata.obs['seurat_clusters'].unique())

        results['seurat'] = {
            'n_clusters': n_seurat_clusters,
            'resolution': seurat_res
        }

        print(f"Seurat: Found {n_seurat_clusters} clusters")
        pandas2ri.deactivate()

    except Exception as e:
        print(f"  Seurat clustering failed: {e}")
        print("  Make sure Seurat is installed: install.packages('Seurat')")
        adata.obs['seurat_clusters'] = 'NA'
        results['seurat'] = {'error': str(e)}

# Celda clustering (R-based)
if run_celda_bool:
    print("Running Celda clustering...")
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.packages import importr

        pandas2ri.activate()

        # Get counts
        if 'counts' in adata.layers:
            counts = adata.layers['counts']
        else:
            counts = adata.X

        # Convert to dense if sparse
        if hasattr(counts, 'toarray'):
            counts_dense = counts.toarray()
        else:
            counts_dense = counts

        # Create R matrix (genes x cells)
        counts_df = pd.DataFrame(
            counts_dense.T,
            index=adata.var_names,
            columns=adata.obs_names
        )

        # Parse celda_L and celda_K
        if celda_L_str.lower() == 'auto':
            celda_L_val = 'NULL'
        else:
            celda_L_val = int(celda_L_str)

        if celda_K_str.lower() == 'auto':
            celda_K_val = 'NULL'
        else:
            celda_K_val = int(celda_K_str)

        ro.r("""
        run_celda_clustering <- function(counts_df, L, K) {
            library(celda)

            # Convert to matrix
            counts_mat <- as.matrix(counts_df)

            # If L and K are NULL, use recursive splitting
            if (is.null(L) || is.null(K)) {
                # Use simple celda_C for cell clustering only
                # Automatically determine number of clusters
                result <- celda_C(counts_mat, K = 5, verbose = FALSE)
                clusters <- celdaClusters(result)$z
            } else {
                # Use celda_CG for both cell and gene clustering
                result <- celda_CG(counts_mat, L = L, K = K, verbose = FALSE)
                clusters <- celdaClusters(result)$z
            }

            return(as.character(clusters))
        }
        """)

        celda_clusters = ro.r['run_celda_clustering'](
            pandas2ri.py2rpy(counts_df),
            ro.r('NULL') if celda_L_str.lower() == 'auto' else celda_L_val,
            ro.r('NULL') if celda_K_str.lower() == 'auto' else celda_K_val
        )

        # Add to adata
        adata.obs['celda_clusters'] = pd.Categorical(list(celda_clusters))
        n_celda_clusters = len(adata.obs['celda_clusters'].unique())

        results['celda'] = {
            'n_clusters': n_celda_clusters,
            'L': celda_L_str,
            'K': celda_K_str
        }

        print(f"Celda: Found {n_celda_clusters} clusters")
        pandas2ri.deactivate()

    except Exception as e:
        print(f"  Celda clustering failed: {e}")
        print("  Make sure celda is installed: BiocManager::install('celda')")
        adata.obs['celda_clusters'] = 'NA'
        results['celda'] = {'error': str(e)}

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
if run_seurat_bool and 'seurat_clusters' in adata.obs:
    cluster_df['seurat'] = adata.obs['seurat_clusters'].values
if run_celda_bool and 'celda_clusters' in adata.obs:
    cluster_df['celda'] = adata.obs['celda_clusters'].values
cluster_df.to_csv('cluster_assignments.csv')

# Create plots
print("Generating clustering plots...")

# Determine number of plots needed
n_cluster_methods = int(run_leiden_bool) + int(run_louvain_bool)
if run_seurat_bool and 'seurat' in results and 'n_clusters' in results.get('seurat', {}):
    n_cluster_methods += 1
if run_celda_bool and 'celda' in results and 'n_clusters' in results.get('celda', {}):
    n_cluster_methods += 1

has_umap = 'X_umap' in adata.obsm
has_pca = 'X_pca' in adata.obsm

# Create figure with appropriate size
n_rows = max(2, n_cluster_methods)
n_cols = 3 if has_umap else 2
fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
if n_rows == 1:
    axes = axes.reshape(1, -1)
if n_cluster_methods == 0:
    n_rows = 2  # Ensure at least 2 rows for layout

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

# Seurat plots
if run_seurat_bool and 'seurat' in results and 'n_clusters' in results.get('seurat', {}):
    # Cluster distribution
    ax = axes[plot_row, 0]
    cluster_counts = adata.obs['seurat_clusters'].value_counts().sort_index()
    ax.bar(range(len(cluster_counts)), cluster_counts.values, color='forestgreen', alpha=0.7)
    ax.set_xlabel('Seurat Cluster')
    ax.set_ylabel('Number of cells')
    ax.set_title(f"Seurat Clustering (n={results['seurat']['n_clusters']})")
    ax.set_xticks(range(len(cluster_counts)))
    ax.set_xticklabels(cluster_counts.index, rotation=45)

    # PCA with clusters
    if has_pca:
        ax = axes[plot_row, 1]
        sc.pl.pca(adata, color='seurat_clusters', ax=ax, show=False, title='Seurat clusters (PCA)', legend_loc='right margin')

    # UMAP with clusters
    if has_umap:
        ax = axes[plot_row, 2] if has_pca else axes[plot_row, 1]
        sc.pl.umap(adata, color='seurat_clusters', ax=ax, show=False, title='Seurat clusters (UMAP)', legend_loc='right margin')

    plot_row += 1

# Celda plots
if run_celda_bool and 'celda' in results and 'n_clusters' in results.get('celda', {}):
    # Cluster distribution
    ax = axes[plot_row, 0]
    cluster_counts = adata.obs['celda_clusters'].value_counts().sort_index()
    ax.bar(range(len(cluster_counts)), cluster_counts.values, color='purple', alpha=0.7)
    ax.set_xlabel('Celda Cluster')
    ax.set_ylabel('Number of cells')
    ax.set_title(f"Celda Clustering (n={results['celda']['n_clusters']})")
    ax.set_xticks(range(len(cluster_counts)))
    ax.set_xticklabels(cluster_counts.index, rotation=45)

    # PCA with clusters
    if has_pca:
        ax = axes[plot_row, 1]
        sc.pl.pca(adata, color='celda_clusters', ax=ax, show=False, title='Celda clusters (PCA)', legend_loc='right margin')

    # UMAP with clusters
    if has_umap:
        ax = axes[plot_row, 2] if has_pca else axes[plot_row, 1]
        sc.pl.umap(adata, color='celda_clusters', ax=ax, show=False, title='Celda clusters (UMAP)', legend_loc='right margin')

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
  Run Seurat: {run_seurat_bool}
  Run Celda: {run_celda_bool}
  Leiden resolution: {leiden_res}
  Louvain resolution: {louvain_res}
  Seurat resolution: {seurat_res}
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

if run_seurat_bool and 'seurat' in results:
    if 'n_clusters' in results['seurat']:
        summary += f"""Seurat SNN Clustering Results:
  Number of clusters: {results['seurat']['n_clusters']}
  Resolution: {results['seurat']['resolution']}

  Cluster sizes:
"""
        for cluster_id in sorted(adata.obs['seurat_clusters'].unique()):
            count = (adata.obs['seurat_clusters'] == cluster_id).sum()
            pct = 100 * count / n_cells
            summary += f"    Cluster {cluster_id}: {count} cells ({pct:.1f}%)\\n"
        summary += "\\n"
    else:
        summary += f"""Seurat Clustering: FAILED
  Error: {results['seurat'].get('error', 'Unknown error')}

"""

if run_celda_bool and 'celda' in results:
    if 'n_clusters' in results['celda']:
        summary += f"""Celda Clustering Results:
  Number of clusters: {results['celda']['n_clusters']}
  L (cell modules): {results['celda']['L']}
  K (gene modules): {results['celda']['K']}

  Cluster sizes:
"""
        for cluster_id in sorted(adata.obs['celda_clusters'].unique()):
            count = (adata.obs['celda_clusters'] == cluster_id).sum()
            pct = 100 * count / n_cells
            summary += f"    Cluster {cluster_id}: {count} cells ({pct:.1f}%)\\n"
        summary += "\\n"
    else:
        summary += f"""Celda Clustering: FAILED
  Error: {results['celda'].get('error', 'Unknown error')}

"""

summary += f"""
Data saved to: clustered.h5ad
Cluster assignments saved to: cluster_assignments.csv

Available cluster annotations:
"""
if run_leiden_bool:
    summary += "  - adata.obs['leiden']\\n"
if run_louvain_bool:
    summary += "  - adata.obs['louvain']\\n"
if run_seurat_bool and 'seurat_clusters' in adata.obs:
    summary += "  - adata.obs['seurat_clusters']\\n"
if run_celda_bool and 'celda_clusters' in adata.obs:
    summary += "  - adata.obs['celda_clusters']\\n"

with open('clustering_summary.txt', 'w') as f:
    f.write(summary)

print(summary)
print("Clustering complete!")
    '''
}
