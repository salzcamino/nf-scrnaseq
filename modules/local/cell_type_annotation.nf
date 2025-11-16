process CELL_TYPE_ANNOTATION {
    tag "annotation"
    label 'process_medium'
    publishDir "${params.outdir}/annotation", mode: params.publish_dir_mode

    input:
    path adata
    val annotation_method
    val celltypist_model
    val marker_file
    val cluster_key

    output:
    path "annotated.h5ad", emit: adata
    path "cell_type_predictions.csv", emit: predictions
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

print("Loading clustered data...")
adata = sc.read_h5ad('!{adata}')

n_cells = adata.n_obs
n_genes = adata.n_vars
print(f"Input: {n_cells} cells, {n_genes} genes")

# Parse parameters
annotation_method = '!{annotation_method}'
celltypist_model_param = '!{celltypist_model}'
marker_file_param = '!{marker_file}'
cluster_key_param = '!{cluster_key}'

# Determine cluster key
if cluster_key_param == 'auto':
    if 'leiden' in adata.obs.columns:
        cluster_key = 'leiden'
    elif 'louvain' in adata.obs.columns:
        cluster_key = 'louvain'
    else:
        cluster_cols = [col for col in adata.obs.columns if 'cluster' in col.lower()]
        cluster_key = cluster_cols[0] if cluster_cols else 'leiden'
else:
    cluster_key = cluster_key_param

print(f"Using cluster key: {cluster_key}")
n_clusters = len(adata.obs[cluster_key].unique())
print(f"Number of clusters: {n_clusters}")

# Initialize results
predictions_df = None
method_used = annotation_method

# Method 1: CellTypist (pre-trained models)
if annotation_method == 'celltypist':
    print("Using CellTypist for cell type annotation...")
    try:
        import celltypist
        from celltypist import models

        # Download model if needed
        if celltypist_model_param == 'Immune_All_Low.pkl':
            print(f"Using model: {celltypist_model_param}")
            models.download_models(model=celltypist_model_param)
            model = models.Model.load(model=celltypist_model_param)
        elif celltypist_model_param == 'Immune_All_High.pkl':
            print(f"Using model: {celltypist_model_param}")
            models.download_models(model=celltypist_model_param)
            model = models.Model.load(model=celltypist_model_param)
        elif Path(celltypist_model_param).exists():
            print(f"Loading custom model: {celltypist_model_param}")
            model = models.Model.load(model=celltypist_model_param)
        else:
            print(f"Downloading model: {celltypist_model_param}")
            models.download_models(model=celltypist_model_param)
            model = models.Model.load(model=celltypist_model_param)

        # CellTypist requires normalized, log-transformed data
        # It expects raw counts in adata.X or .raw
        print("Running CellTypist prediction...")

        # Annotate cells
        predictions = celltypist.annotate(
            adata,
            model=model,
            majority_voting=True
        )

        # Get results
        adata = predictions.to_adata()

        # Extract predictions
        if 'predicted_labels' in adata.obs.columns:
            adata.obs['predicted_cell_type'] = adata.obs['predicted_labels']
        if 'majority_voting' in adata.obs.columns:
            adata.obs['cluster_cell_type'] = adata.obs['majority_voting']

        # Get confidence scores
        if 'conf_score' in adata.obs.columns:
            adata.obs['cell_type_score'] = adata.obs['conf_score']
        else:
            adata.obs['cell_type_score'] = 1.0

        # Create predictions dataframe
        predictions_df = pd.DataFrame({
            'cell': adata.obs_names,
            'predicted_type': adata.obs['predicted_cell_type'].values,
            'confidence': adata.obs['cell_type_score'].values
        })

        print(f"CellTypist annotation complete")
        print(f"  Cell types found: {adata.obs['predicted_cell_type'].nunique()}")

    except ImportError:
        print("WARNING: CellTypist not installed. Falling back to marker scoring.")
        print("Install with: pip install celltypist")
        annotation_method = 'marker_scoring'
        method_used = 'marker_scoring (CellTypist not available)'
    except Exception as e:
        print(f"WARNING: CellTypist failed: {e}")
        print("Falling back to marker scoring method.")
        annotation_method = 'marker_scoring'
        method_used = f'marker_scoring (CellTypist error: {str(e)[:50]})'

# Method 2: Marker gene scoring (fallback or primary)
if annotation_method == 'marker_scoring' or predictions_df is None:
    print("Using marker gene scoring for cell type annotation...")

    # Define default marker genes
    default_markers = {
        'T_cells': ['CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B', 'TRAC'],
        'B_cells': ['CD19', 'MS4A1', 'CD79A', 'CD79B', 'PAX5'],
        'NK_cells': ['NKG7', 'GNLY', 'KLRD1', 'KLRF1', 'NCAM1'],
        'Monocytes': ['CD14', 'LYZ', 'FCGR3A', 'MS4A7', 'CSF1R'],
        'Dendritic_cells': ['FCER1A', 'CLEC10A', 'CD1C', 'ITGAX'],
        'Macrophages': ['CD68', 'CD163', 'MRC1', 'MARCO'],
        'Neutrophils': ['S100A8', 'S100A9', 'FCGR3B', 'CSF3R'],
        'Platelets': ['PPBP', 'PF4', 'GP9', 'ITGA2B'],
        'Erythrocytes': ['HBA1', 'HBA2', 'HBB', 'GYPA'],
        'Epithelial': ['EPCAM', 'KRT8', 'KRT18', 'KRT19', 'CDH1'],
        'Fibroblasts': ['COL1A1', 'COL1A2', 'DCN', 'LUM', 'PDGFRA'],
        'Endothelial': ['PECAM1', 'VWF', 'CDH5', 'KDR', 'FLT1'],
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
            markers = default_markers
    else:
        markers = default_markers

    # Filter markers to genes present in data
    available_genes = set(adata.var_names)
    filtered_markers = {}
    for cell_type, genes in markers.items():
        present_genes = [g for g in genes if g in available_genes]
        if present_genes:
            filtered_markers[cell_type] = present_genes
            print(f"  {cell_type}: {len(present_genes)}/{len(genes)} markers found")

    if len(filtered_markers) == 0:
        print("WARNING: No marker genes found in dataset!")
        adata.obs['predicted_cell_type'] = 'Unknown'
        adata.obs['cell_type_score'] = 0.0
        adata.obs['cluster_cell_type'] = 'Unknown'
        predictions_df = pd.DataFrame({
            'cell': adata.obs_names,
            'predicted_type': ['Unknown'] * n_cells,
            'confidence': [0.0] * n_cells
        })
    else:
        # Score cells for each cell type
        cell_type_scores = {}
        for cell_type, genes in filtered_markers.items():
            if len(genes) >= 2:
                sc.tl.score_genes(adata, genes, score_name=f'{cell_type}_score', use_raw=False)
                cell_type_scores[cell_type] = adata.obs[f'{cell_type}_score'].values
            elif len(genes) == 1:
                gene_idx = list(adata.var_names).index(genes[0])
                if hasattr(adata.X, 'toarray'):
                    scores = np.asarray(adata.X[:, gene_idx].toarray()).flatten()
                else:
                    scores = np.asarray(adata.X[:, gene_idx]).flatten()
                adata.obs[f'{cell_type}_score'] = scores
                cell_type_scores[cell_type] = scores

        # Assign cell types
        if cell_type_scores:
            scores_array = np.array(list(cell_type_scores.values())).T
            cell_type_names = list(cell_type_scores.keys())
            max_scores = np.max(scores_array, axis=1)
            max_indices = np.argmax(scores_array, axis=1)
            predicted_types = [cell_type_names[i] for i in max_indices]

            # Mark low-confidence
            score_threshold = np.percentile(max_scores, 25)
            for i, score in enumerate(max_scores):
                if score < score_threshold:
                    predicted_types[i] = 'Low_confidence'

            adata.obs['predicted_cell_type'] = pd.Categorical(predicted_types)
            adata.obs['cell_type_score'] = max_scores

            predictions_df = pd.DataFrame({
                'cell': adata.obs_names,
                'predicted_type': predicted_types,
                'confidence': max_scores
            })

            # Add cluster-level annotation
            cluster_types = {}
            for cluster in adata.obs[cluster_key].unique():
                cluster_mask = adata.obs[cluster_key] == cluster
                type_counts = adata.obs[cluster_mask]['predicted_cell_type'].value_counts()
                cluster_types[str(cluster)] = type_counts.index[0]

            adata.obs['cluster_cell_type'] = adata.obs[cluster_key].astype(str).map(cluster_types)

# Summarize by cluster
print("Summarizing annotations by cluster...")
cluster_annotations = []

for cluster in sorted(adata.obs[cluster_key].unique(), key=lambda x: str(x)):
    cluster_mask = adata.obs[cluster_key] == cluster
    cluster_cells = adata.obs[cluster_mask]

    type_counts = cluster_cells['predicted_cell_type'].value_counts()
    dominant_type = type_counts.index[0]
    dominant_fraction = type_counts.iloc[0] / len(cluster_cells)
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

# Save predictions
if predictions_df is not None:
    predictions_df.to_csv('cell_type_predictions.csv', index=False)

# Create plots
print("Generating annotation plots...")
fig = plt.figure(figsize=(16, 12))

# Plot 1: UMAP colored by predicted cell type
ax1 = plt.subplot(2, 2, 1)
if 'X_umap' in adata.obsm:
    sc.pl.umap(adata, color='predicted_cell_type', ax=ax1, show=False, legend_loc='right margin',
               legend_fontsize=7, title='Predicted Cell Types (per cell)')
else:
    sc.pl.pca(adata, color='predicted_cell_type', ax=ax1, show=False,
              title='Predicted Cell Types (per cell)')

# Plot 2: UMAP colored by cluster cell type
ax2 = plt.subplot(2, 2, 2)
if 'cluster_cell_type' in adata.obs.columns and 'X_umap' in adata.obsm:
    sc.pl.umap(adata, color='cluster_cell_type', ax=ax2, show=False, legend_loc='right margin',
               legend_fontsize=7, title='Predicted Cell Types (per cluster)')
elif 'X_umap' in adata.obsm:
    sc.pl.umap(adata, color=cluster_key, ax=ax2, show=False, legend_loc='on data',
               title=f'Clusters ({cluster_key})')
else:
    sc.pl.pca(adata, color=cluster_key, ax=ax2, show=False,
              title=f'Clusters ({cluster_key})')

# Plot 3: Cell type distribution
ax3 = plt.subplot(2, 2, 3)
type_counts = adata.obs['predicted_cell_type'].value_counts()
colors = plt.cm.tab20(np.linspace(0, 1, len(type_counts)))
bars = ax3.barh(range(len(type_counts)), type_counts.values, color=colors, alpha=0.8)
ax3.set_yticks(range(len(type_counts)))
ax3.set_yticklabels(type_counts.index, fontsize=9)
ax3.set_xlabel('Number of Cells')
ax3.set_title('Cell Type Distribution')
for i, bar in enumerate(bars):
    width = bar.get_width()
    ax3.text(width, bar.get_y() + bar.get_height()/2,
            f' {int(width)} ({100*width/n_cells:.1f}%)',
            va='center', fontsize=8)

# Plot 4: Confidence score distribution
ax4 = plt.subplot(2, 2, 4)
ax4.hist(adata.obs['cell_type_score'], bins=50, color='steelblue', alpha=0.7, edgecolor='black')
ax4.axvline(adata.obs['cell_type_score'].median(), color='red', linestyle='--',
            label=f'Median: {adata.obs["cell_type_score"].median():.3f}')
ax4.set_xlabel('Confidence Score')
ax4.set_ylabel('Number of Cells')
ax4.set_title('Cell Type Confidence Scores')
ax4.legend()

plt.suptitle(f'Cell Type Annotation Results ({method_used})', fontsize=14, fontweight='bold')
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

Method Used: {method_used}
"""

if method_used.startswith('celltypist'):
    summary += f"CellTypist Model: {celltypist_model_param}\\n"

summary += f"""
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
Mean Confidence Score: {adata.obs['cell_type_score'].mean():.3f}
Median Confidence Score: {adata.obs['cell_type_score'].median():.3f}

Output Files:
  - annotated.h5ad: AnnData with cell type annotations
  - cell_type_predictions.csv: Per-cell type predictions and confidence
  - cluster_annotations.csv: Cluster-level annotations
  - annotation_plots.pdf: Visualization plots

Annotations stored in adata.obs:
  - 'predicted_cell_type': Per-cell type prediction
  - 'cluster_cell_type': Cluster-level type assignment
  - 'cell_type_score': Confidence score
"""

with open('annotation_summary.txt', 'w') as f:
    f.write(summary)

print(summary)
print("Cell type annotation complete!")
    '''
}
