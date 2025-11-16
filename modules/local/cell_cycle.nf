process CELL_CYCLE_SCORING {
    tag "cell_cycle"
    label 'process_low'
    publishDir "${params.outdir}/cell_cycle", mode: params.publish_dir_mode

    input:
    path adata
    val regress_cell_cycle

    output:
    path "cell_cycle_scored.h5ad", emit: adata
    path "cell_cycle_scores.csv", emit: scores
    path "cell_cycle_plots.pdf", emit: plots
    path "cell_cycle_summary.txt", emit: summary

    shell:
    '''
    #!/usr/bin/env python3

    import scanpy as sc
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import warnings
    warnings.filterwarnings('ignore')

    # Parameters
    regress_cc_param = '!{regress_cell_cycle}' == 'true'

    # Load data
    adata = sc.read_h5ad('!{adata}')

    summary_lines = []
    summary_lines.append("Cell Cycle Scoring Summary")
    summary_lines.append("=" * 50)
    summary_lines.append(f"Input cells: {adata.n_obs}")
    summary_lines.append(f"Input genes: {adata.n_vars}")
    summary_lines.append(f"Regress cell cycle: {regress_cc_param}")
    summary_lines.append("")

    # Cell cycle gene sets (from Tirosh et al., 2016)
    # S phase genes
    s_genes = [
        'MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG',
        'GINS2', 'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'MLF1IP',
        'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76',
        'SLBP', 'CCNE2', 'UBR7', 'POLD3', 'MSH2', 'ATAD2', 'RAD51',
        'RRM2', 'CDC45', 'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM',
        'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8'
    ]

    # G2M phase genes
    g2m_genes = [
        'HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A',
        'NDC80', 'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF',
        'TACC3', 'FAM64A', 'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB',
        'BUB1', 'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP',
        'CDCA3', 'HN1', 'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1',
        'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8', 'ECT2', 'KIF23', 'HMMR',
        'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5', 'CENPE', 'CTCF',
        'NEK2', 'G2E3', 'GAS2L3', 'CBX5', 'CENPA'
    ]

    # Check gene availability
    s_genes_available = [g for g in s_genes if g in adata.var_names]
    g2m_genes_available = [g for g in g2m_genes if g in adata.var_names]

    summary_lines.append("Cell Cycle Gene Sets")
    summary_lines.append("-" * 50)
    summary_lines.append(f"S phase genes available: {len(s_genes_available)}/{len(s_genes)}")
    summary_lines.append(f"G2M phase genes available: {len(g2m_genes_available)}/{len(g2m_genes)}")

    if len(s_genes_available) < 5 or len(g2m_genes_available) < 5:
        summary_lines.append("\\nWARNING: Too few cell cycle genes found in dataset")
        summary_lines.append("Cell cycle scoring may not be reliable")
        summary_lines.append("")

    # Store original expression matrix if we're going to regress
    if regress_cc_param:
        adata_original = adata.copy()

    # Score cell cycle phases
    summary_lines.append("\\nScoring Cell Cycle Phases")
    summary_lines.append("-" * 50)

    try:
        sc.tl.score_genes_cell_cycle(
            adata,
            s_genes=s_genes_available,
            g2m_genes=g2m_genes_available
        )
        scoring_success = True
        summary_lines.append("Cell cycle scoring completed successfully")
    except Exception as e:
        summary_lines.append(f"WARNING: Cell cycle scoring failed: {str(e)}")
        summary_lines.append("Attempting manual scoring...")

        # Manual scoring as fallback
        try:
            if len(s_genes_available) >= 2:
                sc.tl.score_genes(adata, s_genes_available, score_name='S_score')
            else:
                adata.obs['S_score'] = 0

            if len(g2m_genes_available) >= 2:
                sc.tl.score_genes(adata, g2m_genes_available, score_name='G2M_score')
            else:
                adata.obs['G2M_score'] = 0

            # Assign phase based on scores
            def assign_phase(row):
                if row['S_score'] > row['G2M_score'] and row['S_score'] > 0.1:
                    return 'S'
                elif row['G2M_score'] > row['S_score'] and row['G2M_score'] > 0.1:
                    return 'G2M'
                else:
                    return 'G1'

            adata.obs['phase'] = adata.obs.apply(assign_phase, axis=1)
            scoring_success = True
            summary_lines.append("Manual scoring completed")
        except Exception as e2:
            summary_lines.append(f"Manual scoring also failed: {str(e2)}")
            scoring_success = False
            adata.obs['S_score'] = 0
            adata.obs['G2M_score'] = 0
            adata.obs['phase'] = 'Unknown'

    # Phase distribution
    if scoring_success and 'phase' in adata.obs.columns:
        phase_counts = adata.obs['phase'].value_counts()
        summary_lines.append("\\nCell Cycle Phase Distribution:")
        for phase, count in phase_counts.items():
            pct = 100 * count / adata.n_obs
            summary_lines.append(f"  {phase}: {count} cells ({pct:.1f}%)")

        # Score statistics
        summary_lines.append("\\nScore Statistics:")
        summary_lines.append(f"  S score: mean={adata.obs['S_score'].mean():.3f}, std={adata.obs['S_score'].std():.3f}")
        summary_lines.append(f"  G2M score: mean={adata.obs['G2M_score'].mean():.3f}, std={adata.obs['G2M_score'].std():.3f}")
    summary_lines.append("")

    # Optional: Regress out cell cycle effects
    if regress_cc_param and scoring_success:
        summary_lines.append("Cell Cycle Regression")
        summary_lines.append("-" * 50)

        try:
            # Store pre-regression data
            adata.uns['pre_cc_regression'] = {
                'S_score_mean': float(adata.obs['S_score'].mean()),
                'G2M_score_mean': float(adata.obs['G2M_score'].mean())
            }

            # Regress out cell cycle scores
            sc.pp.regress_out(adata, ['S_score', 'G2M_score'])

            # Re-scale after regression
            sc.pp.scale(adata, max_value=10)

            summary_lines.append("Cell cycle effects regressed out successfully")
            summary_lines.append("Note: Expression matrix has been modified")
            summary_lines.append("Consider re-running PCA and downstream analyses")

            adata.uns['cell_cycle_regressed'] = True
        except Exception as e:
            summary_lines.append(f"WARNING: Regression failed: {str(e)}")
            summary_lines.append("Continuing without regression")
            adata.uns['cell_cycle_regressed'] = False
    else:
        adata.uns['cell_cycle_regressed'] = False
        if regress_cc_param:
            summary_lines.append("Regression skipped due to scoring issues")
        else:
            summary_lines.append("Cell cycle regression: disabled")

    summary_lines.append("")

    # Store cell cycle info
    adata.uns['cell_cycle'] = {
        's_genes_used': s_genes_available,
        'g2m_genes_used': g2m_genes_available,
        'scoring_successful': scoring_success
    }

    # Save scores
    cc_scores = pd.DataFrame({
        'cell_id': adata.obs_names,
        'S_score': adata.obs['S_score'],
        'G2M_score': adata.obs['G2M_score'],
        'phase': adata.obs['phase']
    })
    cc_scores.to_csv('cell_cycle_scores.csv', index=False)

    # Visualization
    with PdfPages('cell_cycle_plots.pdf') as pdf:
        # Plot 1: Phase distribution
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # Bar plot of phase distribution
        ax = axes[0, 0]
        if 'phase' in adata.obs.columns:
            phase_counts = adata.obs['phase'].value_counts()
            colors = {'G1': '#1f77b4', 'S': '#ff7f0e', 'G2M': '#2ca02c', 'Unknown': '#7f7f7f'}
            phase_colors = [colors.get(p, '#7f7f7f') for p in phase_counts.index]
            phase_counts.plot(kind='bar', ax=ax, color=phase_colors)
            ax.set_title('Cell Cycle Phase Distribution')
            ax.set_xlabel('Phase')
            ax.set_ylabel('Number of Cells')
            ax.tick_params(axis='x', rotation=0)

            # Add percentages on bars
            for i, (phase, count) in enumerate(phase_counts.items()):
                pct = 100 * count / adata.n_obs
                ax.text(i, count + adata.n_obs * 0.01, f'{pct:.1f}%',
                       ha='center', va='bottom')

        # Score distributions
        ax = axes[0, 1]
        ax.hist(adata.obs['S_score'], bins=30, alpha=0.7, label='S score', color='#ff7f0e')
        ax.hist(adata.obs['G2M_score'], bins=30, alpha=0.7, label='G2M score', color='#2ca02c')
        ax.set_title('Cell Cycle Score Distributions')
        ax.set_xlabel('Score')
        ax.set_ylabel('Number of Cells')
        ax.legend()
        ax.axvline(x=0, color='black', linestyle='--', alpha=0.5)

        # S vs G2M scatter
        ax = axes[1, 0]
        if 'phase' in adata.obs.columns:
            for phase in adata.obs['phase'].unique():
                mask = adata.obs['phase'] == phase
                color = colors.get(phase, '#7f7f7f')
                ax.scatter(adata.obs.loc[mask, 'S_score'],
                          adata.obs.loc[mask, 'G2M_score'],
                          alpha=0.5, s=10, label=phase, c=color)
        else:
            ax.scatter(adata.obs['S_score'], adata.obs['G2M_score'],
                      alpha=0.5, s=10)
        ax.set_xlabel('S Score')
        ax.set_ylabel('G2M Score')
        ax.set_title('S vs G2M Phase Scores')
        ax.legend()
        ax.axhline(y=0, color='black', linestyle='--', alpha=0.3)
        ax.axvline(x=0, color='black', linestyle='--', alpha=0.3)

        # UMAP colored by phase (if available)
        ax = axes[1, 1]
        if 'X_umap' in adata.obsm and 'phase' in adata.obs.columns:
            for phase in ['G1', 'S', 'G2M']:
                if phase in adata.obs['phase'].values:
                    mask = adata.obs['phase'] == phase
                    color = colors.get(phase, '#7f7f7f')
                    ax.scatter(adata.obsm['X_umap'][mask, 0],
                              adata.obsm['X_umap'][mask, 1],
                              alpha=0.5, s=10, label=phase, c=color)
            ax.set_xlabel('UMAP1')
            ax.set_ylabel('UMAP2')
            ax.set_title('Cell Cycle Phase on UMAP')
            ax.legend()
        else:
            ax.text(0.5, 0.5, 'UMAP not available',
                   ha='center', va='center', fontsize=12)
            ax.axis('off')

        plt.suptitle('Cell Cycle Analysis', fontsize=14)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()

        # Plot 2: Score heatmap on UMAP (if available)
        if 'X_umap' in adata.obsm:
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))

            # S score on UMAP
            ax = axes[0]
            scatter = ax.scatter(adata.obsm['X_umap'][:, 0],
                                adata.obsm['X_umap'][:, 1],
                                c=adata.obs['S_score'],
                                cmap='YlOrRd', s=10, alpha=0.7)
            plt.colorbar(scatter, ax=ax, label='S Score')
            ax.set_xlabel('UMAP1')
            ax.set_ylabel('UMAP2')
            ax.set_title('S Phase Score')

            # G2M score on UMAP
            ax = axes[1]
            scatter = ax.scatter(adata.obsm['X_umap'][:, 0],
                                adata.obsm['X_umap'][:, 1],
                                c=adata.obs['G2M_score'],
                                cmap='YlGn', s=10, alpha=0.7)
            plt.colorbar(scatter, ax=ax, label='G2M Score')
            ax.set_xlabel('UMAP1')
            ax.set_ylabel('UMAP2')
            ax.set_title('G2M Phase Score')

            plt.suptitle('Cell Cycle Scores on UMAP', fontsize=14)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()

        # Plot 3: Phase by cluster (if clusters available)
        cluster_key = None
        for key in ['leiden', 'louvain', 'seurat_clusters']:
            if key in adata.obs.columns:
                cluster_key = key
                break

        if cluster_key and 'phase' in adata.obs.columns:
            fig, ax = plt.subplots(figsize=(10, 6))

            phase_by_cluster = pd.crosstab(
                adata.obs[cluster_key],
                adata.obs['phase'],
                normalize='index'
            )

            # Reorder columns
            col_order = [c for c in ['G1', 'S', 'G2M'] if c in phase_by_cluster.columns]
            if col_order:
                phase_by_cluster = phase_by_cluster[col_order]

            phase_by_cluster.plot(kind='bar', stacked=True, ax=ax,
                                 color=[colors.get(c, '#7f7f7f') for c in phase_by_cluster.columns])
            ax.set_title(f'Cell Cycle Phase Distribution by {cluster_key.capitalize()} Cluster')
            ax.set_xlabel('Cluster')
            ax.set_ylabel('Proportion')
            ax.legend(title='Phase', bbox_to_anchor=(1.05, 1), loc='upper left')
            ax.tick_params(axis='x', rotation=45)

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()

    # Save data
    adata.write('cell_cycle_scored.h5ad')

    summary_lines.append("Output Files")
    summary_lines.append("-" * 50)
    summary_lines.append("- cell_cycle_scored.h5ad: AnnData with cell cycle scores")
    summary_lines.append("- cell_cycle_scores.csv: Cell cycle scores and phase assignments")
    summary_lines.append("- cell_cycle_plots.pdf: Cell cycle visualization plots")
    summary_lines.append("- cell_cycle_summary.txt: This summary file")

    with open('cell_cycle_summary.txt', 'w') as f:
        f.write("\\n".join(summary_lines))

    print(f"Cell cycle scoring complete. Phases: {adata.obs['phase'].value_counts().to_dict()}")
    '''
}
