process DOUBLET_DECONTAM {
    tag "Doublet detection and decontamination"
    label 'process_medium'

    publishDir "${params.outdir}/doublet_decontam", mode: params.publish_dir_mode

    input:
    path adata
    val run_scrublet
    val run_scdblfinder
    val run_decontx
    val scrublet_threshold
    val expected_doublet_rate

    output:
    path "doublet_scored.h5ad", emit: adata
    path "doublet_scores.csv", emit: scores
    path "doublet_plots.pdf", emit: plots
    path "doublet_summary.txt", emit: summary
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
    import warnings
    warnings.filterwarnings('ignore')

    # Parse parameters
    run_scrublet = '!{run_scrublet}'.lower() == 'true'
    run_scdblfinder = '!{run_scdblfinder}'.lower() == 'true'
    run_decontx = '!{run_decontx}'.lower() == 'true'
    scrublet_threshold = float('!{scrublet_threshold}') if '!{scrublet_threshold}' != 'auto' else None
    expected_doublet_rate = float('!{expected_doublet_rate}')

    # Set plotting parameters
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=80, facecolor='white', frameon=False)
    plt.rcParams['figure.figsize'] = (8, 6)

    # Load data
    print("Loading data...")
    adata = sc.read_h5ad('!{adata}')
    n_cells = adata.n_obs

    print(f"Processing {n_cells} cells...")

    # Initialize results dictionary
    results = {}

    # ========== SCRUBLET DOUBLET DETECTION ==========
    if run_scrublet:
        print("\\nRunning Scrublet doublet detection...")
        try:
            import scrublet as scr

            # Run Scrublet
            scrub = scr.Scrublet(adata.X, expected_doublet_rate=expected_doublet_rate)
            doublet_scores, predicted_doublets = scrub.scrub_doublets(
                min_counts=2,
                min_cells=3,
                min_gene_variability_pct=85,
                n_prin_comps=30
            )

            # Add to AnnData
            adata.obs['scrublet_score'] = doublet_scores
            adata.obs['scrublet_predicted_doublet'] = predicted_doublets

            # Store results
            results['scrublet'] = {
                'n_doublets': int(predicted_doublets.sum()),
                'pct_doublets': float(100 * predicted_doublets.sum() / n_cells),
                'threshold': float(scrub.threshold_) if hasattr(scrub, 'threshold_') else float(scrublet_threshold) if scrublet_threshold else 0.0
            }

            print(f"  Scrublet: {results['scrublet']['n_doublets']} doublets ({results['scrublet']['pct_doublets']:.1f}%)")

        except Exception as e:
            print(f"  Scrublet failed: {e}")
            adata.obs['scrublet_score'] = 0.0
            adata.obs['scrublet_predicted_doublet'] = False
            results['scrublet'] = {'error': str(e)}

    # ========== SCDBLFINDER (R) ==========
    if run_scdblfinder:
        print("\\nRunning scDblFinder...")
        try:
            import anndata2ri
            import rpy2.robjects as ro
            from rpy2.robjects import pandas2ri

            # Activate converters
            anndata2ri.activate()
            pandas2ri.activate()

            # Convert to SingleCellExperiment
            ro.r("""
            library(scDblFinder)
            library(SingleCellExperiment)

            run_scdblfinder <- function(adata) {
                # Convert AnnData to SCE
                sce <- as.SingleCellExperiment(adata)

                # Run scDblFinder
                sce <- scDblFinder(sce)

                # Extract scores
                scores <- data.frame(
                    scdblfinder_score = sce$scDblFinder.score,
                    scdblfinder_class = sce$scDblFinder.class
                )
                return(scores)
            }
            """)

            # Run scDblFinder
            scdblfinder_results = ro.r['run_scdblfinder'](adata)
            scdblfinder_df = pandas2ri.rpy2py(scdblfinder_results)

            # Add to AnnData
            adata.obs['scdblfinder_score'] = scdblfinder_df['scdblfinder_score'].values
            adata.obs['scdblfinder_class'] = scdblfinder_df['scdblfinder_class'].values
            adata.obs['scdblfinder_predicted_doublet'] = (scdblfinder_df['scdblfinder_class'] == 'doublet').values

            n_doublets = adata.obs['scdblfinder_predicted_doublet'].sum()
            results['scdblfinder'] = {
                'n_doublets': int(n_doublets),
                'pct_doublets': float(100 * n_doublets / n_cells)
            }

            print(f"  scDblFinder: {results['scdblfinder']['n_doublets']} doublets ({results['scdblfinder']['pct_doublet']:.1f}%)")

            anndata2ri.deactivate()
            pandas2ri.deactivate()

        except Exception as e:
            print(f"  scDblFinder failed: {e}")
            adata.obs['scdblfinder_score'] = 0.0
            adata.obs['scdblfinder_predicted_doublet'] = False
            results['scdblfinder'] = {'error': str(e)}

    # ========== DECONTX (R) ==========
    if run_decontx:
        print("\\nRunning DecontX for ambient RNA estimation...")
        try:
            import anndata2ri
            import rpy2.robjects as ro
            from rpy2.robjects import pandas2ri

            anndata2ri.activate()
            pandas2ri.activate()

            ro.r("""
            library(celda)
            library(SingleCellExperiment)

            run_decontx <- function(adata) {
                # Convert to SCE
                sce <- as.SingleCellExperiment(adata)

                # Run DecontX
                sce <- decontX(sce)

                # Extract contamination scores
                scores <- data.frame(
                    decontx_contamination = colData(sce)$decontX_contamination
                )
                return(scores)
            }
            """)

            # Run DecontX
            decontx_results = ro.r['run_decontx'](adata)
            decontx_df = pandas2ri.rpy2py(decontx_results)

            # Add to AnnData
            adata.obs['decontx_contamination'] = decontx_df['decontx_contamination'].values

            results['decontx'] = {
                'mean_contamination': float(adata.obs['decontx_contamination'].mean()),
                'median_contamination': float(adata.obs['decontx_contamination'].median())
            }

            print(f"  DecontX: mean contamination = {results['decontx']['mean_contamination']:.2%}")

            anndata2ri.deactivate()
            pandas2ri.deactivate()

        except Exception as e:
            print(f"  DecontX failed: {e}")
            adata.obs['decontx_contamination'] = 0.0
            results['decontx'] = {'error': str(e)}

    # ========== SIMPLE AMBIENT RNA ESTIMATION ==========
    # As backup/complement to DecontX
    print("\\nEstimating ambient RNA (simple method)...")
    try:
        # Estimate based on low-count cells
        low_count_threshold = np.percentile(adata.obs['total_counts'], 10)
        low_count_mask = adata.obs['total_counts'] < low_count_threshold

        if low_count_mask.sum() > 10:
            # Ambient profile from low-count cells
            ambient_profile = np.array(adata[low_count_mask, :].X.mean(axis=0)).flatten()

            # Estimate contamination per cell
            contamination_scores = []
            for i in range(adata.n_obs):
                cell_profile = np.array(adata.X[i, :]).flatten()
                # Correlation with ambient profile
                if cell_profile.sum() > 0:
                    corr = np.corrcoef(cell_profile, ambient_profile)[0, 1]
                    contamination_scores.append(max(0, corr))
                else:
                    contamination_scores.append(0)

            adata.obs['ambient_rna_score'] = contamination_scores
            results['ambient_simple'] = {
                'mean_score': float(np.mean(contamination_scores)),
                'median_score': float(np.median(contamination_scores))
            }
        else:
            adata.obs['ambient_rna_score'] = 0.0
            results['ambient_simple'] = {'note': 'Too few low-count cells for estimation'}

    except Exception as e:
        print(f"  Simple ambient estimation failed: {e}")
        adata.obs['ambient_rna_score'] = 0.0
        results['ambient_simple'] = {'error': str(e)}

    # ========== CREATE VISUALIZATIONS ==========
    print("\\nCreating visualizations...")
    with PdfPages('doublet_plots.pdf') as pdf:

        # Plot 1: Doublet scores from different methods
        n_methods = sum([run_scrublet, run_scdblfinder])
        if n_methods > 0:
            fig, axes = plt.subplots(1, max(n_methods, 1), figsize=(6*n_methods, 5))
            if n_methods == 1:
                axes = [axes]

            plot_idx = 0
            if run_scrublet and 'scrublet_score' in adata.obs:
                axes[plot_idx].hist(adata.obs['scrublet_score'], bins=50, edgecolor='black', alpha=0.7)
                axes[plot_idx].set_xlabel('Scrublet Score')
                axes[plot_idx].set_ylabel('Number of cells')
                axes[plot_idx].set_title('Scrublet Doublet Scores')
                if 'scrublet' in results and 'threshold' in results['scrublet']:
                    axes[plot_idx].axvline(x=results['scrublet']['threshold'], color='red', linestyle='--', label='Threshold')
                    axes[plot_idx].legend()
                plot_idx += 1

            if run_scdblfinder and 'scdblfinder_score' in adata.obs:
                axes[plot_idx].hist(adata.obs['scdblfinder_score'], bins=50, edgecolor='black', alpha=0.7, color='orange')
                axes[plot_idx].set_xlabel('scDblFinder Score')
                axes[plot_idx].set_ylabel('Number of cells')
                axes[plot_idx].set_title('scDblFinder Doublet Scores')
                plot_idx += 1

            plt.tight_layout()
            pdf.savefig(fig)
            plt.close()

        # Plot 2: Contamination scores
        if run_decontx or 'ambient_rna_score' in adata.obs:
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))

            if 'decontx_contamination' in adata.obs:
                axes[0].hist(adata.obs['decontx_contamination'], bins=50, edgecolor='black', alpha=0.7, color='green')
                axes[0].set_xlabel('DecontX Contamination')
                axes[0].set_ylabel('Number of cells')
                axes[0].set_title('DecontX Contamination Scores')

            if 'ambient_rna_score' in adata.obs:
                axes[1].hist(adata.obs['ambient_rna_score'], bins=50, edgecolor='black', alpha=0.7, color='purple')
                axes[1].set_xlabel('Ambient RNA Score')
                axes[1].set_ylabel('Number of cells')
                axes[1].set_title('Ambient RNA Correlation Scores')

            plt.tight_layout()
            pdf.savefig(fig)
            plt.close()

        # Plot 3: Scatter plots - doublet scores vs QC metrics
        if run_scrublet and 'scrublet_score' in adata.obs:
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))

            # Scrublet vs total counts
            scatter = axes[0].scatter(adata.obs['total_counts'],
                                     adata.obs['scrublet_score'],
                                     c=adata.obs['n_genes_by_counts'],
                                     s=5, alpha=0.5, cmap='viridis')
            axes[0].set_xlabel('Total counts')
            axes[0].set_ylabel('Scrublet score')
            axes[0].set_title('Doublet Score vs Total Counts')
            plt.colorbar(scatter, ax=axes[0], label='N genes')

            # Scrublet vs MT percentage
            axes[1].scatter(adata.obs['pct_counts_mt'],
                           adata.obs['scrublet_score'],
                           s=5, alpha=0.5)
            axes[1].set_xlabel('MT %')
            axes[1].set_ylabel('Scrublet score')
            axes[1].set_title('Doublet Score vs MT %')

            plt.tight_layout()
            pdf.savefig(fig)
            plt.close()

    # ========== SAVE RESULTS ==========
    print("\\nSaving results...")

    # Save AnnData with scores
    adata.write('doublet_scored.h5ad', compression='gzip')

    # Save scores table
    score_columns = [col for col in adata.obs.columns if 'doublet' in col.lower() or 'contam' in col.lower() or 'ambient' in col.lower()]
    if score_columns:
        adata.obs[score_columns].to_csv('doublet_scores.csv')
    else:
        pd.DataFrame().to_csv('doublet_scores.csv')

    # Write summary
    with open('doublet_summary.txt', 'w') as f:
        f.write("Doublet Detection and Decontamination Summary\\n")
        f.write("=============================================\\n\\n")
        f.write(f"Total cells analyzed: {n_cells:,}\\n\\n")

        if run_scrublet and 'scrublet' in results:
            f.write("Scrublet Doublet Detection:\\n")
            if 'error' in results['scrublet']:
                f.write(f"  Error: {results['scrublet']['error']}\\n")
            else:
                f.write(f"  Predicted doublets: {results['scrublet']['n_doublets']:,} ({results['scrublet']['pct_doublets']:.2f}%)\\n")
                f.write(f"  Threshold: {results['scrublet']['threshold']:.4f}\\n")
            f.write("\\n")

        if run_scdblfinder and 'scdblfinder' in results:
            f.write("scDblFinder Doublet Detection:\\n")
            if 'error' in results['scdblfinder']:
                f.write(f"  Error: {results['scdblfinder']['error']}\\n")
            else:
                f.write(f"  Predicted doublets: {results['scdblfinder']['n_doublets']:,} ({results['scdblfinder']['pct_doublets']:.2f}%)\\n")
            f.write("\\n")

        if run_decontx and 'decontx' in results:
            f.write("DecontX Contamination Estimation:\\n")
            if 'error' in results['decontx']:
                f.write(f"  Error: {results['decontx']['error']}\\n")
            else:
                f.write(f"  Mean contamination: {results['decontx']['mean_contamination']:.2%}\\n")
                f.write(f"  Median contamination: {results['decontx']['median_contamination']:.2%}\\n")
            f.write("\\n")

        if 'ambient_simple' in results:
            f.write("Simple Ambient RNA Estimation:\\n")
            if 'error' in results['ambient_simple']:
                f.write(f"  Error: {results['ambient_simple']['error']}\\n")
            elif 'note' in results['ambient_simple']:
                f.write(f"  {results['ambient_simple']['note']}\\n")
            else:
                f.write(f"  Mean score: {results['ambient_simple']['mean_score']:.4f}\\n")
                f.write(f"  Median score: {results['ambient_simple']['median_score']:.4f}\\n")
            f.write("\\n")

    # Write versions
    versions = {
        'DOUBLET_DECONTAM': {
            'scanpy': sc.__version__,
            'pandas': pd.__version__,
            'numpy': np.__version__
        }
    }

    try:
        import scrublet
        versions['DOUBLET_DECONTAM']['scrublet'] = scrublet.__version__
    except:
        pass

    try:
        import rpy2
        versions['DOUBLET_DECONTAM']['rpy2'] = rpy2.__version__
    except:
        pass

    with open('versions.yml', 'w') as f:
        yaml.dump(versions, f)

    print("\\nDoublet detection and decontamination completed successfully!")
    '''

    stub:
    '''
    touch doublet_scored.h5ad
    touch doublet_scores.csv
    touch doublet_plots.pdf
    touch doublet_summary.txt
    echo "DOUBLET_DECONTAM:" > versions.yml
    '''
}
