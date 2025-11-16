process CELL_TYPE_ANNOTATION {
    tag "annotation"
    label 'process_medium'
    publishDir "${params.outdir}/annotation", mode: params.publish_dir_mode

    input:
    path adata
    val marker_file
    val cluster_key
    val score_method

    output:
    path "annotated.h5ad", emit: adata
    path "cell_type_scores.csv", emit: scores
    path "cluster_annotations.csv", emit: annotations
    path "annotation_plots.pdf", emit: plots
    path "annotation_summary.txt", emit: summary

    shell:
    '''
#!/usr/bin/env python3
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import json

print("Loading data with DE results...")
adata = sc.read_h5ad('!{adata}')

n_cells = adata.n_obs
n_genes = adata.n_vars
print(f"Input: {n_cells} cells, {n_genes} genes")

# Parse parameters
marker_file_param = '!{marker_file}'
cluster_key_param = '!{cluster_key}'
score_method = '!{score_method}'

# Determine cluster key
if cluster_key_param == 'auto':
    if 'leiden' in adata.obs.columns:
        cluster_key = 'leiden'
    elif 'louvain' in adata.obs.columns:
        cluster_key = 'louvain'
    else:
        cluster_key = [col for col in adata.obs.columns if 'cluster' in col.lower()][0]
else:
    cluster_key = cluster_key_param

print(f"Using cluster key: {cluster_key}")
n_clusters = len(adata.obs[cluster_key].unique())
print(f"Number of clusters: {n_clusters}")

# Define default marker genes for common cell types
# These are canonical markers used in single-cell analysis
default_markers = {
    'T_cells': ['CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B', 'TRAC', 'TRBC1', 'TRBC2'],
    'B_cells': ['CD19', 'MS4A1', 'CD79A', 'CD79B', 'PAX5', 'BANK1', 'BLK'],
    'NK_cells': ['NKG7', 'GNLY', 'KLRD1', 'KLRF1', 'KLRB1', 'NCAM1', 'FCGR3A'],
    'Monocytes': ['CD14', 'LYZ', 'FCGR3A', 'MS4A7', 'ITGAM', 'CSF1R', 'CD68'],
    'Dendritic_cells': ['FCER1A', 'CLEC10A', 'CD1C', 'ITGAX', 'HLA-DRA', 'HLA-DRB1'],
    'Macrophages': ['CD68', 'CD163', 'MRC1', 'MARCO', 'MSR1', 'MERTK'],
    'Neutrophils': ['S100A8', 'S100A9', 'FCGR3B', 'CSF3R', 'CXCR2', 'FPR1'],
    'Platelets': ['PPBP', 'PF4', 'GP9', 'ITGA2B', 'SELP'],
    'Erythrocytes': ['HBA1', 'HBA2', 'HBB', 'GYPA', 'GYPB', 'ALAS2'],
    'Epithelial': ['EPCAM', 'KRT8', 'KRT18', 'KRT19', 'CDH1', 'CLDN4'],
    'Fibroblasts': ['COL1A1', 'COL1A2', 'DCN', 'LUM', 'PDGFRA', 'FAP'],
    'Endothelial': ['PECAM1', 'VWF', 'CDH5', 'KDR', 'FLT1', 'ENG'],
    'Smooth_muscle': ['ACTA2', 'TAGLN', 'MYH11', 'CNN1', 'DES'],
}

# Load custom markers if provided
if marker_file_param != 'default' and Path(marker_file_param).exists():
    print(f"Loading custom markers from {marker_file_param}")
    if marker_file_param.endswith('.json'):
        with open(marker_file_param, 'r') as f:
            markers = json.load(f)
    elif marker_file_param.endswith('.csv'):
        marker_df = pd.read_csv(marker_file_param)
        markers = {}
        for cell_type in marker_df['cell_type'].unique():
            markers[cell_type] = marker_df[marker_df['cell_type'] == cell_type]['gene'].tolist()
    else:
        print(f"  Warning: Unknown marker file format, using defaults")
        markers = default_markers
else:
    print("Using default marker gene sets")
    markers = default_markers

# Filter markers to only include genes present in the data
available_genes = set(adata.var_names)
filtered_markers = {}
for cell_type, genes in markers.items():
    present_genes = [g for g in genes if g in available_genes]
    if len(present_genes) >= 1:
        filtered_markers[cell_type] = present_genes
        print(f"  {cell_type}: {len(present_genes)}/{len(genes)} markers found")
    else:
        print(f"  {cell_type}: No markers found in data, skipping")

if len(filtered_markers) == 0:
    print("WARNING: No marker genes found in the dataset!")
    print("This may indicate the gene names do not match (e.g., using ensembl IDs vs symbols)")
    print("Creating placeholder annotations...")

    # Create placeholder results
    adata.obs['predicted_cell_type'] = 'Unknown'
    adata.obs['cell_type_score'] = 0.0

    scores_df = pd.DataFrame({'cell': adata.obs_names, 'predicted_type': 'Unknown', 'max_score': 0.0})
    scores_df.to_csv('cell_type_scores.csv', index=False)

    cluster_annotations = pd.DataFrame({
        'cluster': sorted(adata.obs[cluster_key].unique(), key=lambda x: str(x)),
        'predicted_type': 'Unknown',
        'mean_score': 0.0,
        'n_cells': [sum(adata.obs[cluster_key] == c) for c in sorted(adata.obs[cluster_key].unique(), key=lambda x: str(x))]
    })
    cluster_annotations.to_csv('cluster_annotations.csv', index=False)

    # Create simple plot
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.text(0.5, 0.5, 'No marker genes found in dataset\\nGene names may not match marker database',
            ha='center', va='center', fontsize=12, transform=ax.transAxes)
    ax.axis('off')
    plt.savefig('annotation_plots.pdf', dpi=150, bbox_inches='tight')
    plt.close()

    summary = f"""Cell Type Annotation Summary
=======================================
WARNING: No marker genes found in dataset

This typically occurs when:
1. Gene names are Ensembl IDs instead of gene symbols
2. Data is from a non-human/mouse organism
3. Custom marker file is needed

Recommendation: Provide a custom marker file with gene names matching your data.
"""
    with open('annotation_summary.txt', 'w') as f:
        f.write(summary)
    adata.write('annotated.h5ad')
    print("Annotation complete (with warnings)")
    exit(0)

# Score cells for each cell type
print("Scoring cells for each cell type...")
cell_type_scores = {}

for cell_type, genes in filtered_markers.items():
    if len(genes) >= 2:
        # Use score_genes for robust scoring
        sc.tl.score_genes(adata, genes, score_name=f'{cell_type}_score', use_raw=False)
        cell_type_scores[cell_type] = adata.obs[f'{cell_type}_score'].values
        print(f"  Scored {cell_type}")
    elif len(genes) == 1:
        # For single gene, use expression directly
        if genes[0] in adata.var_names:
            gene_idx = list(adata.var_names).index(genes[0])
            if hasattr(adata.X, 'toarray'):
                scores = np.asarray(adata.X[:, gene_idx].toarray()).flatten()
            else:
                scores = np.asarray(adata.X[:, gene_idx]).flatten()
            adata.obs[f'{cell_type}_score'] = scores
            cell_type_scores[cell_type] = scores
            print(f"  Scored {cell_type} (single gene)")

# Create scores dataframe
scores_df = pd.DataFrame(cell_type_scores, index=adata.obs_names)
scores_df.to_csv('cell_type_scores.csv')
print(f"Saved scores for {len(cell_type_scores)} cell types")

# Assign cell types based on highest score
print("Assigning cell types...")
if len(cell_type_scores) > 0:
    scores_array = np.array([cell_type_scores[ct] for ct in cell_type_scores.keys()]).T
    cell_type_names = list(cell_type_scores.keys())

    # Get max score and its index for each cell
    max_scores = np.max(scores_array, axis=1)
    max_indices = np.argmax(scores_array, axis=1)

    # Assign cell types
    predicted_types = [cell_type_names[i] for i in max_indices]

    # Mark low-confidence assignments
    score_threshold = np.percentile(max_scores, 25)
    for i, score in enumerate(max_scores):
        if score < score_threshold:
            predicted_types[i] = 'Low_confidence'

    adata.obs['predicted_cell_type'] = pd.Categorical(predicted_types)
    adata.obs['cell_type_score'] = max_scores
else:
    adata.obs['predicted_cell_type'] = 'Unknown'
    adata.obs['cell_type_score'] = 0.0

# Summarize by cluster
print("Summarizing annotations by cluster...")
cluster_annotations = []

for cluster in sorted(adata.obs[cluster_key].unique(), key=lambda x: str(x)):
    cluster_mask = adata.obs[cluster_key] == cluster
    cluster_cells = adata.obs[cluster_mask]

    # Get most common cell type in cluster
    type_counts = cluster_cells['predicted_cell_type'].value_counts()
    dominant_type = type_counts.index[0]
    dominant_fraction = type_counts.iloc[0] / len(cluster_cells)

    # Get mean score for dominant type
    mean_score = cluster_cells['cell_type_score'].mean()

    cluster_annotations.append({
        'cluster': str(cluster),
        'predicted_type': dominant_type,
        'fraction': dominant_fraction,
        'mean_score': mean_score,
        'n_cells': len(cluster_cells)
    })

cluster_df = pd.DataFrame(cluster_annotations)
cluster_df.to_csv('cluster_annotations.csv', index=False)
print(f"Saved cluster annotations")

# Create cluster-level annotations
adata.obs['cluster_cell_type'] = adata.obs[cluster_key].map(
    dict(zip(cluster_df['cluster'], cluster_df['predicted_type']))
)

# Create plots
print("Generating annotation plots...")
fig = plt.figure(figsize=(16, 12))

# Plot 1: UMAP colored by predicted cell type
ax1 = plt.subplot(2, 2, 1)
if 'X_umap' in adata.obsm:
    sc.pl.umap(adata, color='predicted_cell_type', ax=ax1, show=False, legend_loc='on data',
               legend_fontsize=6, title='Predicted Cell Types (per cell)')
else:
    sc.pl.pca(adata, color='predicted_cell_type', ax=ax1, show=False,
              title='Predicted Cell Types (per cell)')

# Plot 2: UMAP colored by cluster cell type
ax2 = plt.subplot(2, 2, 2)
if 'X_umap' in adata.obsm:
    sc.pl.umap(adata, color='cluster_cell_type', ax=ax2, show=False, legend_loc='on data',
               legend_fontsize=6, title='Predicted Cell Types (per cluster)')
else:
    sc.pl.pca(adata, color='cluster_cell_type', ax=ax2, show=False,
              title='Predicted Cell Types (per cluster)')

# Plot 3: Heatmap of cell type scores by cluster
ax3 = plt.subplot(2, 2, 3)
if len(cell_type_scores) > 0:
    # Calculate mean scores per cluster
    cluster_means = []
    for cluster in sorted(adata.obs[cluster_key].unique(), key=lambda x: str(x)):
        cluster_mask = adata.obs[cluster_key] == cluster
        means = {}
        for ct in cell_type_scores.keys():
            means[ct] = adata.obs[cluster_mask][f'{ct}_score'].mean()
        cluster_means.append(means)

    heatmap_df = pd.DataFrame(cluster_means,
                               index=[str(c) for c in sorted(adata.obs[cluster_key].unique(), key=lambda x: str(x))])

    # Normalize by row for better visualization
    heatmap_norm = (heatmap_df.T - heatmap_df.T.min()) / (heatmap_df.T.max() - heatmap_df.T.min() + 1e-10)

    im = ax3.imshow(heatmap_norm.values, aspect='auto', cmap='YlOrRd')
    ax3.set_xticks(range(len(heatmap_df)))
    ax3.set_xticklabels(heatmap_df.index, rotation=45, ha='right')
    ax3.set_yticks(range(len(heatmap_df.columns)))
    ax3.set_yticklabels(heatmap_df.columns, fontsize=8)
    ax3.set_xlabel('Cluster')
    ax3.set_ylabel('Cell Type')
    ax3.set_title('Cell Type Scores by Cluster (normalized)')
    plt.colorbar(im, ax=ax3, label='Normalized Score')
else:
    ax3.text(0.5, 0.5, 'No scores available', ha='center', va='center')
    ax3.axis('off')

# Plot 4: Cell type distribution
ax4 = plt.subplot(2, 2, 4)
type_counts = adata.obs['predicted_cell_type'].value_counts()
bars = ax4.barh(range(len(type_counts)), type_counts.values, color='steelblue', alpha=0.7)
ax4.set_yticks(range(len(type_counts)))
ax4.set_yticklabels(type_counts.index, fontsize=9)
ax4.set_xlabel('Number of Cells')
ax4.set_title('Cell Type Distribution')

# Add value labels
for i, bar in enumerate(bars):
    width = bar.get_width()
    ax4.text(width, bar.get_y() + bar.get_height()/2,
            f' {int(width)} ({100*width/n_cells:.1f}%)',
            va='center', fontsize=8)

plt.suptitle('Cell Type Annotation Results', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('annotation_plots.pdf', dpi=150, bbox_inches='tight')
plt.close()

# Save annotated data
print("Saving annotated data...")
adata.write('annotated.h5ad')

# Write summary
summary = f"""Cell Type Annotation Summary
=======================================
Input cells: {n_cells}
Input genes: {n_genes}
Cluster key: {cluster_key}
Number of clusters: {n_clusters}

Marker Gene Sets:
  Cell types scored: {len(filtered_markers)}
  Scoring method: {score_method}

Cell Type Distribution:
"""

for cell_type, count in type_counts.items():
    pct = 100 * count / n_cells
    summary += f"  {cell_type}: {count} cells ({pct:.1f}%)\\n"

summary += f"""
Cluster Annotations:
"""

for _, row in cluster_df.iterrows():
    summary += f"  Cluster {row['cluster']}: {row['predicted_type']} ({row['fraction']:.1%}, n={row['n_cells']})\\n"

summary += f"""
Output Files:
  - annotated.h5ad: AnnData with cell type annotations
  - cell_type_scores.csv: Cell type scores for each cell
  - cluster_annotations.csv: Cluster-level annotations
  - annotation_plots.pdf: Visualization plots

Annotations stored in adata.obs:
  - 'predicted_cell_type': Per-cell type prediction
  - 'cluster_cell_type': Cluster-level type assignment
  - '*_score': Scores for each cell type
"""

with open('annotation_summary.txt', 'w') as f:
    f.write(summary)

print(summary)
print("Cell type annotation complete!")
    '''
}
