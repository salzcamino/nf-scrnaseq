process DIFF_EXPRESSION {
    tag "diff_expression"
    label 'process_medium'
    publishDir "${params.outdir}/diff_expression", mode: params.publish_dir_mode

    input:
    path adata
    val cluster_key
    val de_method
    val n_genes
    val min_fold_change
    val min_in_group_fraction
    val max_out_group_fraction

    output:
    path "de_results.h5ad", emit: adata
    path "marker_genes.csv", emit: markers
    path "top_markers_per_cluster.csv", emit: top_markers
    path "de_plots.pdf", emit: plots
    path "de_summary.txt", emit: summary

    shell:
    '''
#!/usr/bin/env python3
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

print("Loading clustered data...")
adata = sc.read_h5ad('!{adata}')

n_cells = adata.n_obs
n_genes = adata.n_vars
print(f"Input: {n_cells} cells, {n_genes} genes")

# Parse parameters
cluster_key_param = '!{cluster_key}'
de_method = '!{de_method}'
n_top_genes = int('!{n_genes}')
min_fc = float('!{min_fold_change}')
min_in_group = float('!{min_in_group_fraction}')
max_out_group = float('!{max_out_group_fraction}')

# Determine cluster key to use
if cluster_key_param == 'auto':
    # Try to find a clustering result
    if 'leiden' in adata.obs.columns:
        cluster_key = 'leiden'
    elif 'louvain' in adata.obs.columns:
        cluster_key = 'louvain'
    elif 'seurat_clusters' in adata.obs.columns:
        cluster_key = 'seurat_clusters'
    elif 'celda_clusters' in adata.obs.columns:
        cluster_key = 'celda_clusters'
    else:
        raise ValueError("No clustering results found. Please run clustering first.")
else:
    cluster_key = cluster_key_param
    if cluster_key not in adata.obs.columns:
        raise ValueError(f"Cluster key '{cluster_key}' not found in adata.obs")

print(f"Using cluster key: {cluster_key}")
n_clusters = len(adata.obs[cluster_key].unique())
print(f"Number of clusters: {n_clusters}")

# Ensure we use log-normalized data (not scaled)
# Check if we have raw counts stored
if 'counts' in adata.layers:
    print("Using log-normalized data for DE (from X, raw counts in layers)")
else:
    print("Using current X matrix for DE")

# Run differential expression analysis
print(f"Running differential expression analysis using {de_method} method...")

# Validate method
valid_methods = ['wilcoxon', 't-test', 't-test_overestim_var', 'logreg']
if de_method not in valid_methods:
    print(f"Warning: Invalid method '{de_method}', using 'wilcoxon'")
    de_method = 'wilcoxon'

sc.tl.rank_genes_groups(
    adata,
    groupby=cluster_key,
    method=de_method,
    pts=True,  # Calculate percentage of cells expressing
    key_added='rank_genes_groups'
)

print(f"Differential expression analysis complete")

# Extract results for all clusters
print("Extracting marker genes...")
results_list = []

for cluster in adata.obs[cluster_key].unique():
    cluster_str = str(cluster)

    # Get DE results for this cluster
    names = adata.uns['rank_genes_groups']['names'][cluster_str]
    scores = adata.uns['rank_genes_groups']['scores'][cluster_str]
    pvals = adata.uns['rank_genes_groups']['pvals'][cluster_str]
    pvals_adj = adata.uns['rank_genes_groups']['pvals_adj'][cluster_str]
    logfoldchanges = adata.uns['rank_genes_groups']['logfoldchanges'][cluster_str]

    # Get percentage expressing
    if 'pts' in adata.uns['rank_genes_groups']:
        pts_in = adata.uns['rank_genes_groups']['pts'][cluster_str]
    else:
        pts_in = np.ones(len(names))

    if 'pts_rest' in adata.uns['rank_genes_groups']:
        pts_out = adata.uns['rank_genes_groups']['pts_rest'][cluster_str]
    else:
        pts_out = np.zeros(len(names))

    for i in range(len(names)):
        results_list.append({
            'cluster': cluster_str,
            'gene': names[i],
            'score': scores[i],
            'log2fc': logfoldchanges[i],
            'pval': pvals[i],
            'pval_adj': pvals_adj[i],
            'pct_in_group': pts_in[i] if isinstance(pts_in, np.ndarray) else pts_in,
            'pct_out_group': pts_out[i] if isinstance(pts_out, np.ndarray) else pts_out
        })

results_df = pd.DataFrame(results_list)

# Save all marker genes
results_df.to_csv('marker_genes.csv', index=False)
print(f"Saved {len(results_df)} marker gene results")

# Filter for significant markers
print("Filtering for significant markers...")
filtered_df = results_df[
    (results_df['pval_adj'] < 0.05) &
    (results_df['log2fc'] > np.log2(min_fc)) &
    (results_df['pct_in_group'] > min_in_group) &
    (results_df['pct_out_group'] < max_out_group)
].copy()

# Get top N markers per cluster
top_markers_list = []
for cluster in sorted(adata.obs[cluster_key].unique(), key=lambda x: str(x)):
    cluster_markers = filtered_df[filtered_df['cluster'] == str(cluster)].copy()
    cluster_markers = cluster_markers.sort_values('score', ascending=False)
    top_n = cluster_markers.head(n_top_genes)
    top_markers_list.append(top_n)

if top_markers_list:
    top_markers_df = pd.concat(top_markers_list, ignore_index=True)
else:
    top_markers_df = filtered_df.head(0)

top_markers_df.to_csv('top_markers_per_cluster.csv', index=False)
print(f"Saved top {n_top_genes} markers per cluster ({len(top_markers_df)} total)")

# Create plots
print("Generating differential expression plots...")

# Determine number of plots
fig = plt.figure(figsize=(16, 12))

# Plot 1: Dot plot of top markers
ax1 = plt.subplot(2, 2, 1)
n_show = min(5, n_top_genes)  # Show top 5 genes per cluster
sc.pl.rank_genes_groups_dotplot(
    adata,
    n_genes=n_show,
    groupby=cluster_key,
    ax=ax1,
    show=False,
    dendrogram=False
)
ax1.set_title(f'Top {n_show} Marker Genes per Cluster')

# Plot 2: Heatmap of top markers
ax2 = plt.subplot(2, 2, 2)
try:
    sc.pl.rank_genes_groups_heatmap(
        adata,
        n_genes=n_show,
        groupby=cluster_key,
        show=False,
        ax=ax2,
        show_gene_labels=True,
        dendrogram=False,
        swap_axes=True
    )
    ax2.set_title(f'Heatmap of Top {n_show} Markers')
except Exception as e:
    print(f"  Warning: Could not create heatmap: {e}")
    ax2.text(0.5, 0.5, 'Heatmap not available', ha='center', va='center')
    ax2.axis('off')

# Plot 3: Violin plot of top overall markers
ax3 = plt.subplot(2, 2, 3)
# Get overall top markers (most significant across all clusters)
if len(filtered_df) > 0:
    overall_top = filtered_df.nlargest(min(6, len(filtered_df)), 'score')['gene'].tolist()
    if len(overall_top) > 0:
        sc.pl.violin(
            adata,
            keys=overall_top[:min(4, len(overall_top))],
            groupby=cluster_key,
            ax=ax3,
            show=False,
            rotation=45
        )
        ax3.set_title('Top Marker Gene Expression')
else:
    ax3.text(0.5, 0.5, 'No significant markers found', ha='center', va='center')
    ax3.axis('off')

# Plot 4: Summary statistics
ax4 = plt.subplot(2, 2, 4)
# Count markers per cluster
markers_per_cluster = filtered_df.groupby('cluster').size().reset_index(name='n_markers')
markers_per_cluster = markers_per_cluster.sort_values('cluster', key=lambda x: x.astype(str))

if len(markers_per_cluster) > 0:
    bars = ax4.bar(
        range(len(markers_per_cluster)),
        markers_per_cluster['n_markers'],
        color='steelblue',
        alpha=0.7
    )
    ax4.set_xlabel('Cluster')
    ax4.set_ylabel('Number of Significant Markers')
    ax4.set_title('Significant Marker Genes per Cluster')
    ax4.set_xticks(range(len(markers_per_cluster)))
    ax4.set_xticklabels(markers_per_cluster['cluster'], rotation=45)

    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}', ha='center', va='bottom', fontsize=8)
else:
    ax4.text(0.5, 0.5, 'No significant markers found', ha='center', va='center')
    ax4.axis('off')

plt.suptitle('Differential Expression Analysis Results', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('de_plots.pdf', dpi=150, bbox_inches='tight')
plt.close()

# Save data with DE results
print("Saving data with differential expression results...")
adata.write('de_results.h5ad')

# Write summary
summary = f"""Differential Expression Analysis Summary
=======================================
Input cells: {n_cells}
Input genes: {n_genes}
Cluster key: {cluster_key}
Number of clusters: {n_clusters}

Parameters:
  Method: {de_method}
  N top genes to report: {n_top_genes}
  Min fold change: {min_fc}
  Min fraction in group: {min_in_group}
  Max fraction out group: {max_out_group}

Results:
  Total marker gene tests: {len(results_df)}
  Significant markers (p_adj < 0.05): {len(results_df[results_df['pval_adj'] < 0.05])}
  Filtered markers (all criteria): {len(filtered_df)}

"""

# Add per-cluster summary
summary += "Markers per Cluster:\\n"
for cluster in sorted(adata.obs[cluster_key].unique(), key=lambda x: str(x)):
    cluster_str = str(cluster)
    n_sig = len(filtered_df[filtered_df['cluster'] == cluster_str])
    cluster_top = filtered_df[filtered_df['cluster'] == cluster_str].nlargest(3, 'score')
    top_genes = ', '.join(cluster_top['gene'].tolist()) if len(cluster_top) > 0 else 'None'
    summary += f"  Cluster {cluster_str}: {n_sig} markers\\n"
    summary += f"    Top genes: {top_genes}\\n"

summary += f"""
Output Files:
  - de_results.h5ad: AnnData with DE results stored
  - marker_genes.csv: All marker gene statistics
  - top_markers_per_cluster.csv: Top {n_top_genes} markers per cluster
  - de_plots.pdf: Visualization plots

Accessing results in Python:
  adata.uns['rank_genes_groups']  # Full DE results
  sc.get.rank_genes_groups_df(adata, group='0')  # Get results for cluster 0
"""

with open('de_summary.txt', 'w') as f:
    f.write(summary)

print(summary)
print("Differential expression analysis complete!")
    '''
}
