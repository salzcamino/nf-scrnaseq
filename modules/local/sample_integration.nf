/*
========================================================================================
    Sample Integration Module
========================================================================================
    Integrates multiple samples into a single AnnData object for joint analysis
*/

process SAMPLE_INTEGRATION {
    tag "sample_integration"
    label 'process_medium'
    publishDir "${params.outdir}/integration", mode: params.publish_dir_mode

    input:
    path adata_files  // Multiple h5ad files
    val sample_ids    // Sample identifiers
    val integration_method  // 'concat', 'harmony', 'scanorama', 'scvi'

    output:
    path "integrated.h5ad", emit: adata
    path "integration_summary.txt", emit: summary
    path "integration_plots.pdf", emit: plots

    shell:
    '''
    #!/usr/bin/env python3
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import warnings
    warnings.filterwarnings('ignore')

    # Parse inputs
    adata_files = "!{adata_files}".split()
    sample_ids_str = "!{sample_ids}"
    integration_method = "!{integration_method}"

    # Parse sample IDs
    if sample_ids_str.strip() in ['', 'auto', 'null']:
        # Auto-generate sample IDs from filenames
        sample_ids = []
        for f in adata_files:
            # Extract sample name from filename
            name = f.replace('.h5ad', '').replace('_qc_filtered', '').replace('_raw', '')
            sample_ids.append(name)
    else:
        # Parse provided sample IDs (comma or space separated)
        sample_ids = [s.strip() for s in sample_ids_str.replace(',', ' ').split()]

    # Ensure we have matching sample IDs
    if len(sample_ids) != len(adata_files):
        # Fall back to auto-generated
        sample_ids = [f"sample_{i+1}" for i in range(len(adata_files))]

    summary_lines = []
    summary_lines.append("=" * 60)
    summary_lines.append("Sample Integration Summary")
    summary_lines.append("=" * 60)
    summary_lines.append("")

    # ========================================
    # LOAD AND PREPARE SAMPLES
    # ========================================
    summary_lines.append("Input Samples")
    summary_lines.append("-" * 60)

    adatas = {}
    total_cells = 0
    all_genes = set()

    for adata_file, sample_id in zip(adata_files, sample_ids):
        adata = sc.read_h5ad(adata_file)

        # Add sample information
        adata.obs['sample'] = sample_id
        adata.obs['sample'] = adata.obs['sample'].astype('category')

        # Store original barcode
        if 'original_barcode' not in adata.obs.columns:
            adata.obs['original_barcode'] = adata.obs_names.copy()

        # Make cell barcodes unique by adding sample prefix
        adata.obs_names = [f"{sample_id}_{bc}" for bc in adata.obs_names]

        adatas[sample_id] = adata
        total_cells += adata.n_obs
        all_genes.update(adata.var_names)

        summary_lines.append(f"  {sample_id}: {adata.n_obs:,} cells, {adata.n_vars:,} genes")

    summary_lines.append(f"\\nTotal input cells: {total_cells:,}")
    summary_lines.append(f"Total unique genes across samples: {len(all_genes):,}")
    summary_lines.append("")

    # ========================================
    # CONCATENATE SAMPLES
    # ========================================
    summary_lines.append("Concatenation")
    summary_lines.append("-" * 60)

    # Find common genes (intersection)
    common_genes = set(list(adatas.values())[0].var_names)
    for sample_id, adata in adatas.items():
        common_genes = common_genes.intersection(set(adata.var_names))

    summary_lines.append(f"Common genes across all samples: {len(common_genes):,}")

    # Subset to common genes
    for sample_id in adatas:
        adatas[sample_id] = adatas[sample_id][:, list(common_genes)]

    # Concatenate
    adata_list = [adatas[sid] for sid in sample_ids]
    adata_combined = sc.concat(
        adata_list,
        join='inner',  # Only keep common genes
        label='sample',
        keys=sample_ids,
        index_unique=None  # Already made unique above
    )

    # Ensure sample column is properly set
    if 'sample' not in adata_combined.obs.columns:
        # Reconstruct from batch column if needed
        adata_combined.obs['sample'] = adata_combined.obs.index.str.split('_').str[0]

    adata_combined.obs['sample'] = adata_combined.obs['sample'].astype('category')

    summary_lines.append(f"Combined dataset: {adata_combined.n_obs:,} cells, {adata_combined.n_vars:,} genes")
    summary_lines.append("")

    # Sample composition
    sample_counts = adata_combined.obs['sample'].value_counts()
    summary_lines.append("Sample composition:")
    for sample, count in sample_counts.items():
        pct = (count / adata_combined.n_obs) * 100
        summary_lines.append(f"  {sample}: {count:,} cells ({pct:.1f}%)")

    summary_lines.append("")

    # ========================================
    # INTEGRATION (if requested)
    # ========================================
    integration_applied = False

    if integration_method != 'concat' and len(sample_ids) > 1:
        summary_lines.append(f"Integration Method: {integration_method}")
        summary_lines.append("-" * 60)

        # First, process the data for integration
        # This typically requires: normalization, HVG selection, scaling, PCA
        summary_lines.append("Preparing data for integration...")

        # Normalize if not already done
        if 'log1p' not in adata_combined.uns:
            sc.pp.normalize_total(adata_combined, target_sum=1e4)
            sc.pp.log1p(adata_combined)
            summary_lines.append("  Applied normalization and log1p transformation")

        # Find highly variable genes
        sc.pp.highly_variable_genes(adata_combined, n_top_genes=2000, batch_key='sample')
        n_hvg = adata_combined.var['highly_variable'].sum()
        summary_lines.append(f"  Selected {n_hvg} highly variable genes")

        # Scale data
        sc.pp.scale(adata_combined, max_value=10)

        # Run PCA
        sc.tl.pca(adata_combined, n_comps=50, use_highly_variable=True)
        summary_lines.append("  Computed PCA (50 components)")

        # Apply integration method
        if integration_method == 'harmony':
            try:
                import harmonypy as hm

                # Run Harmony on PCA embeddings
                ho = hm.run_harmony(
                    adata_combined.obsm['X_pca'],
                    adata_combined.obs,
                    'sample',
                    max_iter_harmony=20
                )

                # Store corrected embeddings
                adata_combined.obsm['X_pca_harmony'] = ho.Z_corr.T
                adata_combined.obsm['X_pca_original'] = adata_combined.obsm['X_pca'].copy()
                adata_combined.obsm['X_pca'] = ho.Z_corr.T

                integration_applied = True
                summary_lines.append("  Applied Harmony integration successfully")

            except ImportError:
                summary_lines.append("  WARNING: harmonypy not installed, skipping Harmony integration")
                summary_lines.append("  Install with: pip install harmonypy")
            except Exception as e:
                summary_lines.append(f"  ERROR: Harmony integration failed: {str(e)}")

        elif integration_method == 'scanorama':
            try:
                import scanorama

                # Prepare data for Scanorama
                datasets = []
                genes_lists = []

                for sample_id in sample_ids:
                    mask = adata_combined.obs['sample'] == sample_id
                    sample_data = adata_combined[mask, adata_combined.var['highly_variable']].X
                    if hasattr(sample_data, 'toarray'):
                        sample_data = sample_data.toarray()
                    datasets.append(sample_data)
                    genes_lists.append(list(adata_combined.var_names[adata_combined.var['highly_variable']]))

                # Run Scanorama
                integrated, corrected, genes = scanorama.correct(
                    datasets,
                    genes_lists,
                    return_dimred=True
                )

                # Combine corrected embeddings
                corrected_combined = np.vstack(integrated)

                # Store in obsm
                adata_combined.obsm['X_scanorama'] = corrected_combined

                integration_applied = True
                summary_lines.append(f"  Applied Scanorama integration successfully")
                summary_lines.append(f"  Integrated embedding shape: {corrected_combined.shape}")

            except ImportError:
                summary_lines.append("  WARNING: scanorama not installed, skipping Scanorama integration")
                summary_lines.append("  Install with: pip install scanorama")
            except Exception as e:
                summary_lines.append(f"  ERROR: Scanorama integration failed: {str(e)}")

        elif integration_method == 'bbknn':
            try:
                import bbknn

                # Run BBKNN
                bbknn.bbknn(adata_combined, batch_key='sample', n_pcs=30)

                integration_applied = True
                summary_lines.append("  Applied BBKNN integration successfully")
                summary_lines.append("  Created batch-balanced nearest neighbors graph")

            except ImportError:
                summary_lines.append("  WARNING: bbknn not installed, skipping BBKNN integration")
                summary_lines.append("  Install with: pip install bbknn")
            except Exception as e:
                summary_lines.append(f"  ERROR: BBKNN integration failed: {str(e)}")

        elif integration_method == 'scvi':
            try:
                import scvi

                # Setup anndata for scVI
                scvi.model.SCVI.setup_anndata(
                    adata_combined,
                    layer=None,
                    batch_key='sample'
                )

                # Train scVI model
                model = scvi.model.SCVI(adata_combined, n_layers=2, n_latent=30)
                model.train(max_epochs=100, early_stopping=True)

                # Get latent representation
                adata_combined.obsm['X_scvi'] = model.get_latent_representation()

                integration_applied = True
                summary_lines.append("  Applied scVI integration successfully")
                summary_lines.append(f"  Latent space: {adata_combined.obsm['X_scvi'].shape[1]} dimensions")

            except ImportError:
                summary_lines.append("  WARNING: scvi-tools not installed, skipping scVI integration")
                summary_lines.append("  Install with: pip install scvi-tools")
            except Exception as e:
                summary_lines.append(f"  ERROR: scVI integration failed: {str(e)}")

        else:
            summary_lines.append(f"  WARNING: Unknown integration method '{integration_method}'")

        # Compute neighbors using integrated representation
        if integration_applied:
            summary_lines.append("\\nComputing post-integration embeddings...")

            if integration_method == 'harmony':
                sc.pp.neighbors(adata_combined, use_rep='X_pca', n_neighbors=15)
            elif integration_method == 'scanorama':
                sc.pp.neighbors(adata_combined, use_rep='X_scanorama', n_neighbors=15)
            elif integration_method == 'scvi':
                sc.pp.neighbors(adata_combined, use_rep='X_scvi', n_neighbors=15)
            elif integration_method == 'bbknn':
                pass  # BBKNN already computes neighbors
            else:
                sc.pp.neighbors(adata_combined, n_neighbors=15)

            sc.tl.umap(adata_combined)
            summary_lines.append("  Computed UMAP embedding")

    else:
        if len(sample_ids) == 1:
            summary_lines.append("Single sample - no integration needed")
        else:
            summary_lines.append(f"Integration Method: {integration_method}")
            summary_lines.append("Samples concatenated without integration correction")

    summary_lines.append("")

    # ========================================
    # STORE INTEGRATION METADATA
    # ========================================
    adata_combined.uns['integration'] = {
        'n_samples': len(sample_ids),
        'sample_ids': sample_ids,
        'method': integration_method,
        'integration_applied': integration_applied,
        'n_common_genes': len(common_genes),
        'total_cells': adata_combined.n_obs
    }

    # ========================================
    # VISUALIZATION
    # ========================================
    with PdfPages('integration_plots.pdf') as pdf:
        # Plot 1: Sample composition
        fig, ax = plt.subplots(figsize=(10, 6))
        colors = plt.cm.Set3(np.linspace(0, 1, len(sample_ids)))
        bars = ax.bar(range(len(sample_ids)), sample_counts.values, color=colors, edgecolor='black')
        ax.set_xticks(range(len(sample_ids)))
        ax.set_xticklabels(sample_ids, rotation=45, ha='right')
        ax.set_ylabel('Number of Cells')
        ax.set_title('Sample Composition')

        # Add cell counts on bars
        for bar, count in zip(bars, sample_counts.values):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + total_cells*0.01,
                   f'{count:,}', ha='center', va='bottom', fontsize=9)

        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()

        # Plot 2: UMAP colored by sample (if available)
        if 'X_umap' in adata_combined.obsm:
            fig, ax = plt.subplots(figsize=(12, 10))

            for i, sample_id in enumerate(sample_ids):
                mask = adata_combined.obs['sample'] == sample_id
                ax.scatter(
                    adata_combined.obsm['X_umap'][mask, 0],
                    adata_combined.obsm['X_umap'][mask, 1],
                    c=[colors[i]],
                    label=sample_id,
                    s=5,
                    alpha=0.6
                )

            ax.set_xlabel('UMAP1')
            ax.set_ylabel('UMAP2')
            ax.set_title(f'UMAP colored by Sample (Integration: {integration_method})')
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', markerscale=4)

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()

            # Plot 3: Density plots per sample
            if len(sample_ids) <= 8:
                n_cols = min(4, len(sample_ids))
                n_rows = int(np.ceil(len(sample_ids) / n_cols))
                fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows))
                axes = axes.flatten() if len(sample_ids) > 1 else [axes]

                for i, sample_id in enumerate(sample_ids):
                    ax = axes[i]
                    mask = adata_combined.obs['sample'] == sample_id

                    # Background (all cells)
                    ax.scatter(
                        adata_combined.obsm['X_umap'][:, 0],
                        adata_combined.obsm['X_umap'][:, 1],
                        c='lightgray',
                        s=1,
                        alpha=0.3
                    )

                    # Highlighted sample
                    ax.scatter(
                        adata_combined.obsm['X_umap'][mask, 0],
                        adata_combined.obsm['X_umap'][mask, 1],
                        c=[colors[i]],
                        s=5,
                        alpha=0.8
                    )

                    ax.set_title(f'{sample_id}\\n({mask.sum():,} cells)')
                    ax.set_xlabel('UMAP1')
                    ax.set_ylabel('UMAP2')

                # Hide empty subplots
                for i in range(len(sample_ids), len(axes)):
                    axes[i].axis('off')

                plt.suptitle('Sample Distribution on UMAP', fontsize=14, y=1.02)
                plt.tight_layout()
                pdf.savefig(fig, bbox_inches='tight')
                plt.close()

        # Plot 4: QC metrics by sample (if available)
        qc_metrics = ['n_genes', 'total_counts', 'pct_counts_mt']
        available_metrics = [m for m in qc_metrics if m in adata_combined.obs.columns]

        if available_metrics:
            fig, axes = plt.subplots(1, len(available_metrics), figsize=(5*len(available_metrics), 5))
            if len(available_metrics) == 1:
                axes = [axes]

            for ax, metric in zip(axes, available_metrics):
                data_to_plot = []
                labels = []
                for sample_id in sample_ids:
                    mask = adata_combined.obs['sample'] == sample_id
                    data_to_plot.append(adata_combined.obs.loc[mask, metric].values)
                    labels.append(sample_id)

                parts = ax.violinplot(data_to_plot, positions=range(len(sample_ids)), showmeans=True, showmedians=True)

                # Color violins
                for i, pc in enumerate(parts['bodies']):
                    pc.set_facecolor(colors[i])
                    pc.set_alpha(0.7)

                ax.set_xticks(range(len(sample_ids)))
                ax.set_xticklabels(labels, rotation=45, ha='right')
                ax.set_ylabel(metric.replace('_', ' ').title())
                ax.set_title(f'{metric.replace("_", " ").title()} by Sample')

            plt.suptitle('QC Metrics Distribution by Sample', fontsize=14, y=1.02)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()

    # ========================================
    # SAVE RESULTS
    # ========================================
    adata_combined.write('integrated.h5ad')

    summary_lines.append("Integration Complete")
    summary_lines.append("-" * 60)
    summary_lines.append(f"Final dataset: {adata_combined.n_obs:,} cells × {adata_combined.n_vars:,} genes")
    summary_lines.append(f"Number of samples: {len(sample_ids)}")
    summary_lines.append(f"Integration method: {integration_method}")
    summary_lines.append(f"Integration applied: {integration_applied}")
    summary_lines.append("")
    summary_lines.append("Output files:")
    summary_lines.append("- integrated.h5ad: Combined AnnData object")
    summary_lines.append("- integration_plots.pdf: Visualization of sample composition and mixing")
    summary_lines.append("- integration_summary.txt: This summary file")
    summary_lines.append("")
    summary_lines.append("The 'sample' column in obs contains sample identity for downstream analysis.")

    with open('integration_summary.txt', 'w') as f:
        f.write("\\n".join(summary_lines))

    print(f"Sample integration complete: {adata_combined.n_obs:,} cells from {len(sample_ids)} samples")
    '''
}
