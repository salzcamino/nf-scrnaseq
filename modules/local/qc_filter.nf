process QC_FILTER {
    tag "QC and filtering"
    label 'process_medium'

    publishDir "${params.outdir}/qc", mode: params.publish_dir_mode

    input:
    path adata
    val min_genes
    val min_cells
    val max_genes
    val max_counts
    val max_pct_mt
    val exclude_mt
    val exclude_ribo

    output:
    path "qc_filtered.h5ad", emit: adata
    path "qc_metrics.csv", emit: metrics
    path "qc_plots.pdf", emit: plots
    path "qc_summary.txt", emit: summary
    path "versions.yml", emit: versions

    shell:
    '''
    #!/usr/bin/env python3

    import sys
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.backends.backend_pdf import PdfPages
    import yaml

    # Set plotting parameters
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=80, facecolor='white', frameon=False)
    plt.rcParams['figure.figsize'] = (8, 6)

    # Load data
    print("Loading data...")
    adata = sc.read_h5ad('!{adata}')
    n_cells_raw = adata.n_obs
    n_genes_raw = adata.n_vars

    print(f"Raw data: {n_cells_raw} cells, {n_genes_raw} genes")

    # Calculate QC metrics
    print("Calculating QC metrics...")

    # Identify mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')

    # Identify ribosomal genes
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt', 'ribo'],
        percent_top=None,
        log1p=False,
        inplace=True
    )

    # Save pre-filter metrics
    qc_metrics = adata.obs[[
        'n_genes_by_counts',
        'total_counts',
        'pct_counts_mt',
        'pct_counts_ribo'
    ]].copy()
    qc_metrics['status'] = 'raw'

    # Create QC plots before filtering
    print("Creating QC plots...")
    with PdfPages('qc_plots.pdf') as pdf:
        # Plot 1: Distribution of counts, genes, and MT percentage
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # Total counts
        axes[0, 0].hist(adata.obs['total_counts'], bins=50, edgecolor='black')
        axes[0, 0].set_xlabel('Total counts')
        axes[0, 0].set_ylabel('Number of cells')
        axes[0, 0].set_title('Total counts per cell')
        axes[0, 0].axvline(x=!{max_counts if max_counts != 'null' else 'np.inf'},
                          color='red', linestyle='--', label='Max threshold')
        axes[0, 0].legend()

        # Genes per cell
        axes[0, 1].hist(adata.obs['n_genes_by_counts'], bins=50, edgecolor='black')
        axes[0, 1].set_xlabel('Number of genes')
        axes[0, 1].set_ylabel('Number of cells')
        axes[0, 1].set_title('Genes per cell')
        axes[0, 1].axvline(x=!{min_genes}, color='red', linestyle='--', label='Min threshold')
        axes[0, 1].axvline(x=!{max_genes}, color='red', linestyle='--', label='Max threshold')
        axes[0, 1].legend()

        # MT percentage
        axes[1, 0].hist(adata.obs['pct_counts_mt'], bins=50, edgecolor='black')
        axes[1, 0].set_xlabel('Mitochondrial %')
        axes[1, 0].set_ylabel('Number of cells')
        axes[1, 0].set_title('Mitochondrial percentage')
        axes[1, 0].axvline(x=!{max_pct_mt}, color='red', linestyle='--', label='Max threshold')
        axes[1, 0].legend()

        # Scatter: counts vs genes
        axes[1, 1].scatter(adata.obs['total_counts'],
                          adata.obs['n_genes_by_counts'],
                          c=adata.obs['pct_counts_mt'],
                          s=1, alpha=0.5, cmap='viridis')
        axes[1, 1].set_xlabel('Total counts')
        axes[1, 1].set_ylabel('Number of genes')
        axes[1, 1].set_title('Counts vs Genes (colored by MT %)')
        plt.colorbar(axes[1, 1].collections[0], ax=axes[1, 1], label='MT %')

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

        # Plot 2: Violin plots
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))

        sc.pl.violin(adata, ['n_genes_by_counts'], ax=axes[0], show=False)
        axes[0].set_title('Genes per cell')

        sc.pl.violin(adata, ['total_counts'], ax=axes[1], show=False)
        axes[1].set_title('Total counts')

        sc.pl.violin(adata, ['pct_counts_mt'], ax=axes[2], show=False)
        axes[2].set_title('MT percentage')

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

        # Plot 3: Scatter plots
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=axes[0], show=False)
        axes[0].set_title('Total counts vs MT %')

        sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax=axes[1], show=False)
        axes[1].set_title('Total counts vs Genes')

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

    # Apply filters
    print("Applying filters...")

    # Cell filtering
    print(f"Filtering cells with <{!{min_genes}} or >{!{max_genes}} genes...")
    sc.pp.filter_cells(adata, min_genes=!{min_genes})
    adata = adata[adata.obs['n_genes_by_counts'] < !{max_genes}, :].copy()

    # Max counts filter
    max_counts_val = !{max_counts if max_counts != 'null' else 'np.inf'}
    if max_counts_val != np.inf:
        print(f"Filtering cells with >{max_counts_val} counts...")
        adata = adata[adata.obs['total_counts'] < max_counts_val, :].copy()

    # MT percentage filter
    print(f"Filtering cells with >{!{max_pct_mt}}% MT genes...")
    adata = adata[adata.obs['pct_counts_mt'] < !{max_pct_mt}, :].copy()

    # Gene filtering
    print(f"Filtering genes expressed in <{!{min_cells}} cells...")
    sc.pp.filter_genes(adata, min_cells=!{min_cells})

    # Exclude MT genes if requested
    if !{exclude_mt}:
        print("Removing mitochondrial genes...")
        adata = adata[:, ~adata.var['mt']].copy()

    # Exclude ribosomal genes if requested
    if !{exclude_ribo}:
        print("Removing ribosomal genes...")
        adata = adata[:, ~adata.var['ribo']].copy()

    n_cells_filtered = adata.n_obs
    n_genes_filtered = adata.n_vars

    print(f"Filtered data: {n_cells_filtered} cells, {n_genes_filtered} genes")
    print(f"Removed {n_cells_raw - n_cells_filtered} cells ({100*(n_cells_raw - n_cells_filtered)/n_cells_raw:.1f}%)")
    print(f"Removed {n_genes_raw - n_genes_filtered} genes ({100*(n_genes_raw - n_genes_filtered)/n_genes_raw:.1f}%)")

    # Save filtered data
    adata.write('qc_filtered.h5ad', compression='gzip')

    # Save QC metrics
    qc_metrics.to_csv('qc_metrics.csv')

    # Write summary
    with open('qc_summary.txt', 'w') as f:
        f.write("QC Filtering Summary\\n")
        f.write("===================\\n\\n")
        f.write(f"Raw data:\\n")
        f.write(f"  Cells: {n_cells_raw:,}\\n")
        f.write(f"  Genes: {n_genes_raw:,}\\n\\n")
        f.write(f"Filtering parameters:\\n")
        f.write(f"  Min genes per cell: {!{min_genes}}\\n")
        f.write(f"  Max genes per cell: {!{max_genes}}\\n")
        f.write(f"  Min cells per gene: {!{min_cells}}\\n")
        if max_counts_val != np.inf:
            f.write(f"  Max counts per cell: {max_counts_val}\\n")
        f.write(f"  Max MT percentage: {!{max_pct_mt}}%\\n")
        f.write(f"  Exclude MT genes: {!{exclude_mt}}\\n")
        f.write(f"  Exclude ribo genes: {!{exclude_ribo}}\\n\\n")
        f.write(f"Filtered data:\\n")
        f.write(f"  Cells: {n_cells_filtered:,}\\n")
        f.write(f"  Genes: {n_genes_filtered:,}\\n\\n")
        f.write(f"Removed:\\n")
        f.write(f"  Cells: {n_cells_raw - n_cells_filtered:,} ({100*(n_cells_raw - n_cells_filtered)/n_cells_raw:.1f}%)\\n")
        f.write(f"  Genes: {n_genes_raw - n_genes_filtered:,} ({100*(n_genes_raw - n_genes_filtered)/n_genes_raw:.1f}%)\\n")

    # Write versions
    versions = {
        'QC_FILTER': {
            'scanpy': sc.__version__,
            'pandas': pd.__version__,
            'numpy': np.__version__,
            'matplotlib': plt.matplotlib.__version__,
            'seaborn': sns.__version__
        }
    }
    with open('versions.yml', 'w') as f:
        yaml.dump(versions, f)

    print("QC filtering completed successfully!")
    '''

    stub:
    '''
    touch qc_filtered.h5ad
    touch qc_metrics.csv
    touch qc_plots.pdf
    touch qc_summary.txt
    echo "QC_FILTER:" > versions.yml
    '''
}
