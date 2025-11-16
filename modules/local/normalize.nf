process NORMALIZE {
    tag "normalize"
    label 'process_medium'
    publishDir "${params.outdir}/normalize", mode: params.publish_dir_mode

    input:
    path adata
    val target_sum
    val log_transform

    output:
    path "normalized.h5ad", emit: adata
    path "normalization_summary.txt", emit: summary

    shell:
    '''
    #!/usr/bin/env python3
    import scanpy as sc
    import numpy as np
    from pathlib import Path

    print("Loading QC-filtered data...")
    adata = sc.read_h5ad('!{adata}')

    n_cells = adata.n_obs
    n_genes = adata.n_vars
    print(f"Input: {n_cells} cells, {n_genes} genes")

    # Store raw counts in a layer
    adata.layers['counts'] = adata.X.copy()

    # Parse target_sum
    target_sum_str = '!{target_sum}'
    if target_sum_str.lower() in ['auto', 'median', 'null']:
        target_sum_val = None  # Scanpy will use median of counts
        target_sum_desc = "median of counts"
    else:
        target_sum_val = float(target_sum_str)
        target_sum_desc = f"{target_sum_val:.0f}"

    # Parse log_transform
    do_log = '!{log_transform}'.lower() == 'true'

    # Normalize to target sum (library size normalization)
    print(f"Normalizing to target sum: {target_sum_desc}")
    sc.pp.normalize_total(adata, target_sum=target_sum_val)

    # Calculate median counts after normalization
    median_counts = np.median(adata.X.sum(axis=1))

    # Log transform
    if do_log:
        print("Log-transforming data (log1p)")
        sc.pp.log1p(adata)
        transform_desc = "log1p"
    else:
        transform_desc = "none"

    # Save normalized data
    print("Saving normalized data...")
    adata.write('normalized.h5ad')

    # Write summary
    summary = f"""Normalization Summary
====================
Input cells: {n_cells}
Input genes: {n_genes}

Normalization:
  Target sum: {target_sum_desc}
  Log transform: {transform_desc}

Result:
  Median library size: {median_counts:.2f}
  Raw counts stored in: adata.layers['counts']
"""

    with open('normalization_summary.txt', 'w') as f:
        f.write(summary)

    print(summary)
    print("Normalization complete!")
    '''
}
