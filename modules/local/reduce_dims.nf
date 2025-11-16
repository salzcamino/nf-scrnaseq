process REDUCE_DIMS {
    tag "reduce_dims"
    label 'process_medium'
    publishDir "${params.outdir}/dim_reduction", mode: params.publish_dir_mode

    input:
    path adata
    val n_pcs
    val n_neighbors
    val run_umap
    val run_tsne
    val umap_min_dist
    val tsne_perplexity

    output:
    path "reduced_dims.h5ad", emit: adata
    path "dim_reduction_plots.pdf", emit: plots
    path "dim_reduction_summary.txt", emit: summary

    shell:
    '''
    #!/usr/bin/env python3
    import scanpy as sc
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    from pathlib import Path

    print("Loading data with HVG annotations...")
    adata = sc.read_h5ad('!{adata}')

    n_cells = adata.n_obs
    n_genes = adata.n_vars
    n_hvg = adata.var['highly_variable'].sum()
    print(f"Input: {n_cells} cells, {n_genes} genes ({n_hvg} HVGs)")

    # Parse parameters
    n_pcs_val = int('!{n_pcs}')
    n_neighbors_val = int('!{n_neighbors}')
    run_umap_bool = '!{run_umap}'.lower() == 'true'
    run_tsne_bool = '!{run_tsne}'.lower() == 'true'
    umap_min_dist_val = float('!{umap_min_dist}')
    tsne_perplexity_val = float('!{tsne_perplexity}')

    # Adjust parameters for small datasets
    if n_cells < n_pcs_val:
        n_pcs_val = n_cells - 1
        print(f"Adjusted n_pcs to {n_pcs_val} (less than n_cells)")

    if n_neighbors_val >= n_cells:
        n_neighbors_val = max(2, n_cells - 1)
        print(f"Adjusted n_neighbors to {n_neighbors_val}")

    if tsne_perplexity_val >= n_cells / 3:
        tsne_perplexity_val = max(5, n_cells / 5)
        print(f"Adjusted t-SNE perplexity to {tsne_perplexity_val}")

    results = {}

    # Scale data (zero center and unit variance)
    print("Scaling data to zero mean and unit variance...")
    sc.pp.scale(adata, max_value=10)

    # PCA
    print(f"Running PCA with {n_pcs_val} components...")
    sc.tl.pca(adata, n_comps=n_pcs_val, svd_solver='arpack')

    # Calculate variance explained
    var_explained = adata.uns['pca']['variance_ratio']
    cumulative_var = np.cumsum(var_explained)

    results['pca'] = {
        'n_components': n_pcs_val,
        'variance_explained_total': float(cumulative_var[-1]),
        'variance_pc1': float(var_explained[0]),
        'variance_pc2': float(var_explained[1]) if len(var_explained) > 1 else 0
    }

    print(f"PCA complete: {results['pca']['variance_explained_total']:.1%} variance explained")

    # Compute neighbors graph (needed for UMAP)
    print(f"Computing neighbor graph with {n_neighbors_val} neighbors...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors_val, n_pcs=min(50, n_pcs_val))

    # UMAP
    if run_umap_bool:
        print(f"Running UMAP (min_dist={umap_min_dist_val})...")
        sc.tl.umap(adata, min_dist=umap_min_dist_val)
        results['umap'] = {'min_dist': umap_min_dist_val}
        print("UMAP complete")

    # t-SNE
    if run_tsne_bool:
        print(f"Running t-SNE (perplexity={tsne_perplexity_val})...")
        sc.tl.tsne(adata, perplexity=tsne_perplexity_val, n_pcs=min(50, n_pcs_val))
        results['tsne'] = {'perplexity': tsne_perplexity_val}
        print("t-SNE complete")

    # Create plots
    print("Generating dimensionality reduction plots...")

    # Determine number of plots
    n_plots = 2  # PCA variance + PC1 vs PC2
    if run_umap_bool:
        n_plots += 1
    if run_tsne_bool:
        n_plots += 1

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    plot_idx = 0

    # Plot 1: PCA variance explained
    ax = axes[plot_idx]
    x = np.arange(1, len(var_explained) + 1)
    ax.bar(x, var_explained * 100, alpha=0.7, label='Individual')
    ax.plot(x, cumulative_var * 100, 'r-', linewidth=2, label='Cumulative')
    ax.set_xlabel('Principal Component')
    ax.set_ylabel('Variance Explained (%)')
    ax.set_title('PCA Variance Explained')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plot_idx += 1

    # Plot 2: PC1 vs PC2
    ax = axes[plot_idx]
    sc.pl.pca(adata, ax=ax, show=False, title='PCA (PC1 vs PC2)')
    plot_idx += 1

    # Plot 3: UMAP (if enabled)
    if run_umap_bool:
        ax = axes[plot_idx]
        sc.pl.umap(adata, ax=ax, show=False, title='UMAP')
        plot_idx += 1

    # Plot 4: t-SNE (if enabled)
    if run_tsne_bool:
        ax = axes[plot_idx]
        sc.pl.tsne(adata, ax=ax, show=False, title='t-SNE')
        plot_idx += 1

    # Hide unused axes
    for i in range(plot_idx, len(axes)):
        axes[i].axis('off')

    plt.suptitle('Dimensionality Reduction Results', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('dim_reduction_plots.pdf', dpi=150, bbox_inches='tight')
    plt.close()

    # Save data
    print("Saving data with dimensionality reductions...")
    adata.write('reduced_dims.h5ad')

    # Write summary
    summary = f"""Dimensionality Reduction Summary
================================
Input cells: {n_cells}
Input genes: {n_genes}
Highly variable genes: {n_hvg}

Parameters:
  N PCs: {n_pcs_val}
  N neighbors: {n_neighbors_val}
  Run UMAP: {run_umap_bool}
  Run t-SNE: {run_tsne_bool}
  UMAP min_dist: {umap_min_dist_val}
  t-SNE perplexity: {tsne_perplexity_val}

PCA Results:
  Components computed: {n_pcs_val}
  Total variance explained: {results['pca']['variance_explained_total']:.1%}
  PC1 variance: {results['pca']['variance_pc1']:.1%}
  PC2 variance: {results['pca']['variance_pc2']:.1%}
"""

    if run_umap_bool:
        summary += f"""
UMAP Results:
  Embedding computed successfully
  Min distance: {umap_min_dist_val}
"""

    if run_tsne_bool:
        summary += f"""
t-SNE Results:
  Embedding computed successfully
  Perplexity: {tsne_perplexity_val}
"""

    summary += f"""
Data saved to: reduced_dims.h5ad
Available embeddings:
  - adata.obsm['X_pca']
"""

    if run_umap_bool:
        summary += "  - adata.obsm['X_umap']\\n"
    if run_tsne_bool:
        summary += "  - adata.obsm['X_tsne']\\n"

    with open('dim_reduction_summary.txt', 'w') as f:
        f.write(summary)

    print(summary)
    print("Dimensionality reduction complete!")
    '''
}
