/*
========================================================================================
    HTML Report Generation Module
========================================================================================
    Generates comprehensive HTML reports consolidating all pipeline outputs
*/

process HTML_REPORT {
    tag "html_report"
    label 'process_low'
    publishDir "${params.outdir}/report", mode: params.publish_dir_mode

    input:
    path adata

    output:
    path "pipeline_report.html", emit: report
    path "report_data/", emit: data_dir
    path "report_plots/", emit: plots_dir

    shell:
    '''
    #!/usr/bin/env python3
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    import base64
    import json
    from pathlib import Path
    from datetime import datetime
    import io
    import warnings
    warnings.filterwarnings('ignore')

    # Create output directories
    Path("report_data").mkdir(exist_ok=True)
    Path("report_plots").mkdir(exist_ok=True)

    # Load AnnData
    adata = sc.read_h5ad("!{adata}")

    # Helper function to save plot and encode as base64
    def save_and_encode_plot(fig, filename):
        """Save plot to file and return base64 encoding."""
        filepath = f"report_plots/{filename}"
        fig.savefig(filepath, dpi=120, bbox_inches='tight', facecolor='white')
        plt.close(fig)

        with open(filepath, "rb") as f:
            return base64.b64encode(f.read()).decode('utf-8')

    # Collect pipeline statistics
    stats = {
        'n_cells': adata.n_obs,
        'n_genes': adata.n_vars,
        'report_time': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
    }

    # Get clustering info
    cluster_keys = [k for k in adata.obs.columns if k in ['leiden', 'louvain', 'seurat_clusters', 'celda_clusters']]
    if cluster_keys:
        stats['primary_clustering'] = cluster_keys[0]
        stats['n_clusters'] = adata.obs[cluster_keys[0]].nunique()
    else:
        stats['n_clusters'] = 'N/A'

    # Get cell type info
    if 'cell_type' in adata.obs.columns:
        stats['n_cell_types'] = adata.obs['cell_type'].nunique()
        cell_type_counts = adata.obs['cell_type'].value_counts().to_dict()
    else:
        stats['n_cell_types'] = 'N/A'
        cell_type_counts = {}

    # Get cell cycle info
    if 'phase' in adata.obs.columns:
        phase_counts = adata.obs['phase'].value_counts().to_dict()
    else:
        phase_counts = {}

    # Generate plots and encode them
    encoded_plots = {}

    # QC Plots
    if 'n_genes' in adata.obs and 'total_counts' in adata.obs:
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))

        # Genes per cell
        axes[0].hist(adata.obs['n_genes'], bins=50, color='#3498db', edgecolor='black', alpha=0.7)
        axes[0].set_xlabel('Number of Genes')
        axes[0].set_ylabel('Number of Cells')
        axes[0].set_title('Genes per Cell')
        axes[0].axvline(adata.obs['n_genes'].median(), color='red', linestyle='--', label=f"Median: {adata.obs['n_genes'].median():.0f}")
        axes[0].legend()

        # Counts per cell
        axes[1].hist(adata.obs['total_counts'], bins=50, color='#27ae60', edgecolor='black', alpha=0.7)
        axes[1].set_xlabel('Total Counts')
        axes[1].set_ylabel('Number of Cells')
        axes[1].set_title('Counts per Cell')
        axes[1].axvline(adata.obs['total_counts'].median(), color='red', linestyle='--', label=f"Median: {adata.obs['total_counts'].median():.0f}")
        axes[1].legend()

        # MT percentage
        if 'pct_counts_mt' in adata.obs:
            axes[2].hist(adata.obs['pct_counts_mt'], bins=50, color='#e74c3c', edgecolor='black', alpha=0.7)
            axes[2].set_xlabel('MT %')
            axes[2].set_ylabel('Number of Cells')
            axes[2].set_title('Mitochondrial Percentage')
            axes[2].axvline(adata.obs['pct_counts_mt'].median(), color='blue', linestyle='--', label=f"Median: {adata.obs['pct_counts_mt'].median():.1f}%")
            axes[2].legend()
        else:
            axes[2].text(0.5, 0.5, 'MT % not available', ha='center', va='center')
            axes[2].set_axis_off()

        plt.tight_layout()
        encoded_plots['qc_distributions'] = save_and_encode_plot(fig, 'qc_distributions.png')

    # UMAP/Clustering plots
    if 'X_umap' in adata.obsm and cluster_keys:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        primary_key = cluster_keys[0]
        umap_coords = adata.obsm['X_umap']

        # Colored by cluster
        scatter = axes[0].scatter(umap_coords[:, 0], umap_coords[:, 1],
                                  c=adata.obs[primary_key].astype('category').cat.codes,
                                  cmap='tab20', s=5, alpha=0.7)
        axes[0].set_xlabel('UMAP1')
        axes[0].set_ylabel('UMAP2')
        axes[0].set_title(f'UMAP colored by {primary_key.title()}')

        # Add cluster labels
        for cluster in adata.obs[primary_key].unique():
            mask = adata.obs[primary_key] == cluster
            center = umap_coords[mask].mean(axis=0)
            axes[0].text(center[0], center[1], str(cluster), fontsize=10, fontweight='bold',
                        ha='center', va='center', bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

        # Colored by cell type if available
        if 'cell_type' in adata.obs.columns:
            cell_types = adata.obs['cell_type'].astype('category')
            scatter2 = axes[1].scatter(umap_coords[:, 0], umap_coords[:, 1],
                                       c=cell_types.cat.codes, cmap='tab20', s=5, alpha=0.7)
            axes[1].set_xlabel('UMAP1')
            axes[1].set_ylabel('UMAP2')
            axes[1].set_title('UMAP colored by Cell Type')

            # Add legend
            for i, ct in enumerate(cell_types.cat.categories[:10]):  # Show top 10
                axes[1].scatter([], [], c=[plt.cm.tab20(i/20)], label=ct, s=50)
            axes[1].legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
        else:
            axes[1].set_axis_off()

        plt.tight_layout()
        encoded_plots['umap_clusters'] = save_and_encode_plot(fig, 'umap_clusters.png')

    # Cell cycle plot
    if phase_counts:
        fig, ax = plt.subplots(figsize=(8, 6))
        phases = ['G1', 'S', 'G2M']
        counts = [phase_counts.get(p, 0) for p in phases]
        colors = ['#3498db', '#27ae60', '#e74c3c']

        bars = ax.bar(phases, counts, color=colors, edgecolor='black', alpha=0.8)
        ax.set_ylabel('Number of Cells')
        ax.set_title('Cell Cycle Phase Distribution')

        # Add percentage labels
        total = sum(counts)
        for bar, count in zip(bars, counts):
            pct = (count / total) * 100
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + total*0.01,
                   f'{pct:.1f}%', ha='center', va='bottom', fontweight='bold')

        plt.tight_layout()
        encoded_plots['cell_cycle'] = save_and_encode_plot(fig, 'cell_cycle.png')

    # Pseudotime plot
    if 'dpt_pseudotime' in adata.obs.columns and 'X_umap' in adata.obsm:
        fig, ax = plt.subplots(figsize=(10, 8))
        scatter = ax.scatter(adata.obsm['X_umap'][:, 0], adata.obsm['X_umap'][:, 1],
                            c=adata.obs['dpt_pseudotime'], cmap='viridis', s=10, alpha=0.8)
        plt.colorbar(scatter, ax=ax, label='Pseudotime')
        ax.set_xlabel('UMAP1')
        ax.set_ylabel('UMAP2')
        ax.set_title('Diffusion Pseudotime on UMAP')
        plt.tight_layout()
        encoded_plots['pseudotime'] = save_and_encode_plot(fig, 'pseudotime.png')

    # Generate HTML report
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>nf-scrnaseq Pipeline Report</title>
    <style>
        :root {{
            --primary-color: #2c3e50;
            --secondary-color: #3498db;
            --bg-color: #ecf0f1;
            --card-bg: #ffffff;
            --text-color: #2c3e50;
            --border-color: #bdc3c7;
        }}

        * {{ margin: 0; padding: 0; box-sizing: border-box; }}

        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background-color: var(--bg-color);
            color: var(--text-color);
            line-height: 1.6;
        }}

        .header {{
            background: linear-gradient(135deg, var(--primary-color), var(--secondary-color));
            color: white;
            padding: 2rem;
            text-align: center;
        }}

        .header h1 {{ font-size: 2.5rem; margin-bottom: 0.5rem; }}
        .header .subtitle {{ font-size: 1.1rem; opacity: 0.9; }}

        .container {{ max-width: 1400px; margin: 0 auto; padding: 2rem; }}

        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 1.5rem;
            margin-bottom: 2rem;
        }}

        .stat-card {{
            background: var(--card-bg);
            padding: 1.5rem;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.05);
            text-align: center;
            transition: transform 0.2s;
        }}

        .stat-card:hover {{ transform: translateY(-5px); }}
        .stat-value {{ font-size: 2.5rem; font-weight: bold; color: var(--secondary-color); }}
        .stat-label {{ font-size: 0.9rem; color: #7f8c8d; text-transform: uppercase; letter-spacing: 1px; }}

        .section {{
            background: var(--card-bg);
            padding: 2rem;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.05);
            margin-bottom: 2rem;
        }}

        .section h2 {{
            color: var(--primary-color);
            margin-bottom: 1.5rem;
            padding-bottom: 0.5rem;
            border-bottom: 3px solid var(--secondary-color);
        }}

        .plot-container {{ text-align: center; margin: 1.5rem 0; }}
        .plot-container img {{ max-width: 100%; height: auto; border-radius: 5px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}

        .table-container {{ overflow-x: auto; margin: 1rem 0; }}
        table {{ width: 100%; border-collapse: collapse; font-size: 0.9rem; }}
        th, td {{ padding: 0.75rem; text-align: left; border-bottom: 1px solid var(--border-color); }}
        th {{ background: var(--primary-color); color: white; font-weight: 600; }}
        tr:nth-child(even) {{ background: #f8f9fa; }}
        tr:hover {{ background: #e8f4f8; }}

        .collapsible {{
            cursor: pointer;
            padding: 1rem;
            background: #f8f9fa;
            border: none;
            text-align: left;
            width: 100%;
            font-size: 1rem;
            font-weight: 500;
            border-radius: 5px;
            margin-bottom: 0.5rem;
        }}

        .collapsible:hover {{ background: #e8f4f8; }}
        .collapsible:after {{ content: '\\25BC'; float: right; }}
        .collapsible.collapsed:after {{ content: '\\25B6'; }}
        .collapsible-content {{ padding: 1rem; background: #fafafa; border-radius: 0 0 5px 5px; margin-bottom: 1rem; }}
        .hidden {{ display: none; }}
        .status-available {{ color: #27ae60; font-weight: bold; }}
        .status-unavailable {{ color: #95a5a6; }}

        .footer {{ text-align: center; padding: 2rem; color: #7f8c8d; font-size: 0.9rem; }}

        @media (max-width: 768px) {{
            .header h1 {{ font-size: 1.8rem; }}
            .container {{ padding: 1rem; }}
            .section {{ padding: 1rem; }}
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>nf-scrnaseq Pipeline Report</h1>
        <div class="subtitle">Single-cell RNA-seq Analysis Results</div>
        <div class="subtitle" style="margin-top: 0.5rem; opacity: 0.8;">Generated: {stats['report_time']}</div>
    </div>

    <div class="container">
        <div class="summary-grid">
            <div class="stat-card">
                <div class="stat-value">{stats['n_cells']:,}</div>
                <div class="stat-label">Total Cells</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{stats['n_genes']:,}</div>
                <div class="stat-label">Total Genes</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{stats['n_clusters']}</div>
                <div class="stat-label">Clusters</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{stats['n_cell_types']}</div>
                <div class="stat-label">Cell Types</div>
            </div>
        </div>

        <div class="section">
            <h2>Quality Control</h2>
            <div class="table-container">
                <table>
                    <thead><tr><th>Metric</th><th>Value</th></tr></thead>
                    <tbody>
                        <tr><td>Cells after QC</td><td>{adata.n_obs:,}</td></tr>
                        <tr><td>Genes after QC</td><td>{adata.n_vars:,}</td></tr>
"""

    # Add QC metrics
    if 'n_genes' in adata.obs:
        html_content += f"                        <tr><td>Mean genes/cell</td><td>{adata.obs['n_genes'].mean():.1f}</td></tr>\\n"
        html_content += f"                        <tr><td>Median genes/cell</td><td>{adata.obs['n_genes'].median():.1f}</td></tr>\\n"
    if 'total_counts' in adata.obs:
        html_content += f"                        <tr><td>Mean counts/cell</td><td>{adata.obs['total_counts'].mean():.1f}</td></tr>\\n"
        html_content += f"                        <tr><td>Median counts/cell</td><td>{adata.obs['total_counts'].median():.1f}</td></tr>\\n"
    if 'pct_counts_mt' in adata.obs:
        html_content += f"                        <tr><td>Mean MT %</td><td>{adata.obs['pct_counts_mt'].mean():.2f}%</td></tr>\\n"
    if 'doublet_score' in adata.obs:
        html_content += f"                        <tr><td>Doublet detection</td><td>Performed</td></tr>\\n"
        if 'predicted_doublet' in adata.obs:
            n_doublets = adata.obs['predicted_doublet'].sum()
            html_content += f"                        <tr><td>Predicted doublets</td><td>{n_doublets} ({n_doublets/adata.n_obs*100:.1f}%)</td></tr>\\n"

    html_content += """
                    </tbody>
                </table>
            </div>
"""

    # Add QC plot
    if 'qc_distributions' in encoded_plots:
        html_content += f"""
            <div class="plot-container">
                <img src="data:image/png;base64,{encoded_plots['qc_distributions']}" alt="QC Distributions">
            </div>
"""

    # Clustering section
    html_content += """
        </div>

        <div class="section">
            <h2>Clustering Results</h2>
"""

    if cluster_keys:
        primary_key = cluster_keys[0]
        cluster_counts = adata.obs[primary_key].value_counts().sort_index()
        html_content += f"""
            <h3>Cluster Distribution ({primary_key.title()})</h3>
            <div class="table-container">
                <table>
                    <thead><tr><th>Cluster</th><th>Number of Cells</th><th>Percentage</th></tr></thead>
                    <tbody>
"""
        for cluster, count in cluster_counts.items():
            pct = (count / adata.n_obs) * 100
            html_content += f"                        <tr><td>{cluster}</td><td>{count:,}</td><td>{pct:.1f}%</td></tr>\\n"
        html_content += """
                    </tbody>
                </table>
            </div>
"""
        if 'umap_clusters' in encoded_plots:
            html_content += f"""
            <div class="plot-container">
                <img src="data:image/png;base64,{encoded_plots['umap_clusters']}" alt="UMAP Clustering">
            </div>
"""
    else:
        html_content += """
            <p class="status-unavailable">Clustering not performed.</p>
"""

    # Cell type section
    html_content += """
        </div>

        <div class="section">
            <h2>Cell Type Annotation</h2>
"""

    if cell_type_counts:
        html_content += """
            <div class="table-container">
                <table>
                    <thead><tr><th>Cell Type</th><th>Number of Cells</th><th>Percentage</th></tr></thead>
                    <tbody>
"""
        for cell_type, count in sorted(cell_type_counts.items(), key=lambda x: -x[1]):
            pct = (count / adata.n_obs) * 100
            html_content += f"                        <tr><td>{cell_type}</td><td>{count:,}</td><td>{pct:.1f}%</td></tr>\\n"
        html_content += """
                    </tbody>
                </table>
            </div>
"""
    else:
        html_content += """
            <p class="status-unavailable">Cell type annotation not performed.</p>
"""

    # Differential Expression section
    html_content += """
        </div>

        <div class="section">
            <h2>Differential Expression</h2>
"""

    if 'rank_genes_groups' in adata.uns:
        html_content += """
            <p><span class="status-available">✓</span> Differential expression analysis completed</p>
"""
        try:
            result = adata.uns['rank_genes_groups']
            groups = result['names'].dtype.names
            n_genes_show = min(10, len(result['names']))

            for group in list(groups)[:6]:
                html_content += f"""
            <button class="collapsible collapsed">Cluster {group} - Top {n_genes_show} Markers</button>
            <div class="collapsible-content hidden">
                <div class="table-container">
                    <table>
                        <thead><tr><th>Rank</th><th>Gene</th><th>Score</th><th>Log2FC</th><th>Adj. P-value</th></tr></thead>
                        <tbody>
"""
                for i in range(n_genes_show):
                    gene = result['names'][i][group]
                    score = result['scores'][i][group]
                    logfc = result['logfoldchanges'][i][group] if 'logfoldchanges' in result else 'N/A'
                    pval = result['pvals_adj'][i][group] if 'pvals_adj' in result else result['pvals'][i][group]
                    logfc_str = f"{logfc:.3f}" if isinstance(logfc, (float, np.floating)) else str(logfc)
                    html_content += f"                            <tr><td>{i+1}</td><td><strong>{gene}</strong></td><td>{score:.3f}</td><td>{logfc_str}</td><td>{pval:.2e}</td></tr>\\n"

                html_content += """
                        </tbody>
                    </table>
                </div>
            </div>
"""
        except Exception as e:
            html_content += f"            <p class='status-unavailable'>Error loading DE results: {str(e)}</p>\\n"
    else:
        html_content += """
            <p class="status-unavailable">Differential expression analysis not performed.</p>
"""

    # Cell cycle section
    html_content += """
        </div>

        <div class="section">
            <h2>Cell Cycle Analysis</h2>
"""

    if phase_counts:
        html_content += """
            <p><span class="status-available">✓</span> Cell cycle scoring completed</p>
            <div class="table-container">
                <table>
                    <thead><tr><th>Phase</th><th>Number of Cells</th><th>Percentage</th></tr></thead>
                    <tbody>
"""
        for phase in ['G1', 'S', 'G2M']:
            if phase in phase_counts:
                count = phase_counts[phase]
                pct = (count / adata.n_obs) * 100
                html_content += f"                        <tr><td>{phase}</td><td>{count:,}</td><td>{pct:.1f}%</td></tr>\\n"
        html_content += """
                    </tbody>
                </table>
            </div>
"""
        if 'cell_cycle' in encoded_plots:
            html_content += f"""
            <div class="plot-container">
                <img src="data:image/png;base64,{encoded_plots['cell_cycle']}" alt="Cell Cycle Distribution">
            </div>
"""
    else:
        html_content += """
            <p class="status-unavailable">Cell cycle analysis not performed.</p>
"""

    # Trajectory section
    html_content += """
        </div>

        <div class="section">
            <h2>Trajectory Analysis</h2>
"""

    if 'dpt_pseudotime' in adata.obs.columns:
        html_content += f"""
            <p><span class="status-available">✓</span> Trajectory analysis completed</p>
            <p><strong>Pseudotime range:</strong> {adata.obs['dpt_pseudotime'].min():.4f} - {adata.obs['dpt_pseudotime'].max():.4f}</p>
"""
        if 'pseudotime' in encoded_plots:
            html_content += f"""
            <div class="plot-container">
                <img src="data:image/png;base64,{encoded_plots['pseudotime']}" alt="Pseudotime">
            </div>
"""
    else:
        html_content += """
            <p class="status-unavailable">Trajectory analysis not performed.</p>
"""

    # Cell communication section
    html_content += """
        </div>

        <div class="section">
            <h2>Cell-Cell Communication</h2>
"""

    if 'communication_results' in adata.uns:
        html_content += """
            <p><span class="status-available">✓</span> Cell-cell communication analysis completed</p>
            <p>Results are available in the communication output directory.</p>
"""
    else:
        html_content += """
            <p class="status-unavailable">Cell-cell communication analysis not performed or results not stored in AnnData.</p>
"""

    # GSEA section
    html_content += """
        </div>

        <div class="section">
            <h2>Gene Set Enrichment Analysis</h2>
"""

    if 'gsea_results' in adata.uns or any('enrichment' in str(k) for k in adata.uns.keys()):
        html_content += """
            <p><span class="status-available">✓</span> Gene set enrichment analysis completed</p>
            <p>Results are available in the GSEA output directory.</p>
"""
    else:
        html_content += """
            <p class="status-unavailable">Gene set enrichment analysis not performed.</p>
"""

    # Dataset info section
    html_content += """
        </div>

        <div class="section">
            <h2>Dataset Information</h2>
            <button class="collapsible collapsed">AnnData Structure</button>
            <div class="collapsible-content hidden">
                <div class="table-container">
                    <table>
                        <thead><tr><th>Component</th><th>Contents</th></tr></thead>
                        <tbody>
                            <tr><td>obs (cell metadata)</td><td>""" + ", ".join(list(adata.obs.columns)[:15]) + ("..." if len(adata.obs.columns) > 15 else "") + """</td></tr>
                            <tr><td>var (gene metadata)</td><td>""" + ", ".join(list(adata.var.columns)[:10]) + ("..." if len(adata.var.columns) > 10 else "") + """</td></tr>
                            <tr><td>obsm (embeddings)</td><td>""" + ", ".join(list(adata.obsm.keys())) + """</td></tr>
                            <tr><td>uns (unstructured)</td><td>""" + ", ".join(list(adata.uns.keys())[:10]) + ("..." if len(adata.uns.keys()) > 10 else "") + """</td></tr>
                        </tbody>
                    </table>
                </div>
            </div>

            <button class="collapsible collapsed">Available Embeddings</button>
            <div class="collapsible-content hidden">
                <div class="table-container">
                    <table>
                        <thead><tr><th>Embedding</th><th>Dimensions</th></tr></thead>
                        <tbody>
"""

    for key in adata.obsm.keys():
        shape = adata.obsm[key].shape
        html_content += f"                            <tr><td>{key}</td><td>{shape[1]} components</td></tr>\\n"

    html_content += """
                        </tbody>
                    </table>
                </div>
            </div>
        </div>
    </div>

    <div class="footer">
        <p>Generated by <strong>nf-scrnaseq</strong> | Single-cell RNA-seq Analysis Pipeline</p>
        <p>Report generated on """ + stats['report_time'] + """</p>
    </div>

    <script>
        document.querySelectorAll('.collapsible').forEach(button => {
            button.addEventListener('click', function() {
                this.classList.toggle('collapsed');
                this.nextElementSibling.classList.toggle('hidden');
            });
        });
    </script>
</body>
</html>
"""

    # Write HTML file
    with open("pipeline_report.html", 'w') as f:
        f.write(html_content)

    # Save summary data
    summary_data = {
        'statistics': stats,
        'cell_type_counts': cell_type_counts,
        'phase_counts': phase_counts,
        'obs_columns': list(adata.obs.columns),
        'var_columns': list(adata.var.columns),
        'obsm_keys': list(adata.obsm.keys()),
        'uns_keys': list(adata.uns.keys())[:50]
    }

    with open("report_data/summary.json", 'w') as f:
        json.dump(summary_data, f, indent=2, default=str)

    adata.obs.to_csv("report_data/cell_metadata.csv")
    adata.var.to_csv("report_data/gene_metadata.csv")

    print(f"HTML report generated successfully!")
    print(f"Report file: pipeline_report.html")
    print(f"Data directory: report_data/")
    print(f"Plots directory: report_plots/")
    '''
}
