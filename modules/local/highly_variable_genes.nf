process HIGHLY_VARIABLE_GENES {
    tag "hvg"
    label 'process_medium'
    publishDir "${params.outdir}/hvg", mode: params.publish_dir_mode

    input:
    path adata
    val n_top_genes
    val min_mean
    val max_mean
    val min_disp
    val flavor

    output:
    path "hvg_selected.h5ad", emit: adata
    path "hvg_genes.csv", emit: genes
    path "hvg_plots.pdf", emit: plots
    path "hvg_summary.txt", emit: summary

    shell:
    '''
    #!/usr/bin/env python3
    import scanpy as sc
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    from pathlib import Path

    print("Loading normalized data...")
    adata = sc.read_h5ad('!{adata}')

    n_cells = adata.n_obs
    n_genes = adata.n_vars
    print(f"Input: {n_cells} cells, {n_genes} genes")

    # Parse parameters
    n_top_genes_str = '!{n_top_genes}'
    if n_top_genes_str.lower() in ['none', 'null', 'auto']:
        n_top_genes_val = None
    else:
        n_top_genes_val = int(n_top_genes_str)

    min_mean_val = float('!{min_mean}')
    max_mean_val = float('!{max_mean}')
    min_disp_val = float('!{min_disp}')
    flavor_str = '!{flavor}'

    print(f"Identifying highly variable genes...")
    print(f"  Flavor: {flavor_str}")
    if n_top_genes_val:
        print(f"  Top N genes: {n_top_genes_val}")
    else:
        print(f"  Min mean: {min_mean_val}")
        print(f"  Max mean: {max_mean_val}")
        print(f"  Min dispersion: {min_disp_val}")

    # Identify highly variable genes
    if n_top_genes_val:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes_val,
            flavor=flavor_str
        )
    else:
        sc.pp.highly_variable_genes(
            adata,
            min_mean=min_mean_val,
            max_mean=max_mean_val,
            min_disp=min_disp_val,
            flavor=flavor_str
        )

    n_hvg = adata.var['highly_variable'].sum()
    print(f"Found {n_hvg} highly variable genes")

    # Save HVG information
    hvg_df = adata.var[['highly_variable', 'means', 'dispersions', 'dispersions_norm']].copy()
    hvg_df = hvg_df.sort_values('dispersions_norm', ascending=False)
    hvg_df.to_csv('hvg_genes.csv')

    # Create plots
    print("Generating HVG plots...")
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Plot 1: Mean vs dispersion
    ax = axes[0, 0]
    means = adata.var['means']
    dispersions = adata.var['dispersions_norm']
    hvg_mask = adata.var['highly_variable']

    ax.scatter(means[~hvg_mask], dispersions[~hvg_mask], s=3, alpha=0.3, c='gray', label='Non-HVG')
    ax.scatter(means[hvg_mask], dispersions[hvg_mask], s=5, alpha=0.6, c='red', label='HVG')
    ax.set_xlabel('Mean expression')
    ax.set_ylabel('Normalized dispersion')
    ax.set_title(f'Highly Variable Genes (n={n_hvg})')
    ax.legend()
    ax.set_xscale('log')

    # Plot 2: Distribution of means for HVGs
    ax = axes[0, 1]
    ax.hist(means[hvg_mask], bins=50, color='red', alpha=0.7, edgecolor='black')
    ax.set_xlabel('Mean expression')
    ax.set_ylabel('Number of genes')
    ax.set_title('Distribution of mean expression (HVGs)')
    ax.set_xscale('log')

    # Plot 3: Distribution of dispersions for HVGs
    ax = axes[1, 0]
    ax.hist(dispersions[hvg_mask], bins=50, color='blue', alpha=0.7, edgecolor='black')
    ax.set_xlabel('Normalized dispersion')
    ax.set_ylabel('Number of genes')
    ax.set_title('Distribution of normalized dispersion (HVGs)')

    # Plot 4: Top 20 HVGs by dispersion
    ax = axes[1, 1]
    top_hvgs = hvg_df[hvg_df['highly_variable']].head(20)
    y_pos = np.arange(len(top_hvgs))
    ax.barh(y_pos, top_hvgs['dispersions_norm'].values, color='green', alpha=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top_hvgs.index, fontsize=8)
    ax.set_xlabel('Normalized dispersion')
    ax.set_title('Top 20 Highly Variable Genes')
    ax.invert_yaxis()

    plt.suptitle('Highly Variable Gene Selection', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('hvg_plots.pdf', dpi=150, bbox_inches='tight')
    plt.close()

    # Save data with HVG annotations
    print("Saving data with HVG annotations...")
    adata.write('hvg_selected.h5ad')

    # Write summary
    summary = f"""Highly Variable Gene Selection Summary
======================================
Input cells: {n_cells}
Input genes: {n_genes}

Selection parameters:
  Flavor: {flavor_str}
  N top genes: {n_top_genes_val if n_top_genes_val else 'Not used'}
  Min mean: {min_mean_val}
  Max mean: {max_mean_val}
  Min dispersion: {min_disp_val}

Results:
  Highly variable genes: {n_hvg}
  Percentage of total: {100 * n_hvg / n_genes:.1f}%

Top 10 HVGs by normalized dispersion:
"""

    for i, (gene, row) in enumerate(top_hvgs.head(10).iterrows()):
        summary += f"  {i+1}. {gene}: {row['dispersions_norm']:.3f}\\n"

    with open('hvg_summary.txt', 'w') as f:
        f.write(summary)

    print(summary)
    print("HVG selection complete!")
    '''
}
