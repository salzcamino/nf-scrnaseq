process GENE_SET_ENRICHMENT {
    tag "gsea"
    label 'process_medium'
    publishDir "${params.outdir}/gsea", mode: params.publish_dir_mode

    input:
    path adata
    path marker_genes
    val gene_sets
    val n_top_genes
    val cluster_key

    output:
    path "gsea_results.h5ad", emit: adata
    path "enrichment_results.csv", emit: enrichment
    path "pathway_scores.csv", emit: scores
    path "gsea_plots.pdf", emit: plots
    path "gsea_summary.txt", emit: summary

    shell:
    '''
#!/usr/bin/env python3
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

print("Loading annotated data...")
adata = sc.read_h5ad('!{adata}')

n_cells = adata.n_obs
n_genes = adata.n_vars
print(f"Input: {n_cells} cells, {n_genes} genes")

# Parse parameters
gene_sets_param = '!{gene_sets}'
n_top_genes_param = int('!{n_top_genes}')
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
clusters = sorted(adata.obs[cluster_key].unique(), key=lambda x: str(x))
n_clusters = len(clusters)
print(f"Number of clusters: {n_clusters}")

# Load marker genes from DE analysis
print("Loading marker genes...")
markers_df = pd.read_csv('!{marker_genes}')
print(f"Loaded {len(markers_df)} marker gene results")

# Define common gene sets for pathway scoring
# These are canonical pathways from MSigDB Hallmark collection
hallmark_sets = {
    'HALLMARK_INFLAMMATORY_RESPONSE': [
        'IL1B', 'IL6', 'TNF', 'CCL2', 'CXCL8', 'IL1A', 'CCL3', 'CCL4', 'ICAM1', 'NFKBIA'
    ],
    'HALLMARK_INTERFERON_ALPHA_RESPONSE': [
        'ISG15', 'IFI27', 'IFIT1', 'IFIT3', 'MX1', 'OAS1', 'OAS2', 'STAT1', 'IRF7', 'IFITM1'
    ],
    'HALLMARK_INTERFERON_GAMMA_RESPONSE': [
        'STAT1', 'IRF1', 'GBP1', 'CXCL10', 'CXCL9', 'IDO1', 'CIITA', 'TAP1', 'PSMB9', 'HLA-DRA'
    ],
    'HALLMARK_TNFA_SIGNALING_VIA_NFKB': [
        'NFKBIA', 'TNF', 'JUNB', 'ATF3', 'DUSP1', 'ZFP36', 'FOSB', 'JUN', 'RELB', 'BCL2A1'
    ],
    'HALLMARK_APOPTOSIS': [
        'BCL2', 'BAX', 'CASP3', 'CASP8', 'FAS', 'FASLG', 'BID', 'APAF1', 'CYCS', 'DIABLO'
    ],
    'HALLMARK_P53_PATHWAY': [
        'TP53', 'CDKN1A', 'BAX', 'MDM2', 'GADD45A', 'PMAIP1', 'BBC3', 'SESN1', 'DDB2', 'RRM2B'
    ],
    'HALLMARK_OXIDATIVE_PHOSPHORYLATION': [
        'ATP5F1A', 'ATP5F1B', 'COX4I1', 'COX5A', 'NDUFS1', 'NDUFS2', 'SDHA', 'SDHB', 'UQCRC1', 'UQCRC2'
    ],
    'HALLMARK_MYC_TARGETS_V1': [
        'MYC', 'NCL', 'NPM1', 'NME1', 'PA2G4', 'LDHA', 'ENO1', 'PKM', 'HNRNPA1', 'DDX18'
    ],
    'HALLMARK_G2M_CHECKPOINT': [
        'CDC20', 'CCNB1', 'CCNB2', 'BUB1', 'PLK1', 'AURKA', 'AURKB', 'CDK1', 'TOP2A', 'MKI67'
    ],
    'HALLMARK_E2F_TARGETS': [
        'E2F1', 'MCM2', 'MCM3', 'MCM4', 'MCM5', 'MCM6', 'MCM7', 'PCNA', 'RFC4', 'POLE'
    ],
}

# Cell type specific gene sets
cell_type_sets = {
    'T_CELL_ACTIVATION': ['CD69', 'CD25', 'IL2RA', 'ICOS', 'CD40LG', 'CTLA4', 'PDCD1', 'LAG3'],
    'B_CELL_ACTIVATION': ['CD69', 'CD86', 'AICDA', 'BCL6', 'IRF4', 'PRDM1', 'XBP1'],
    'CYTOTOXICITY': ['PRF1', 'GZMA', 'GZMB', 'GZMH', 'GZMK', 'GNLY', 'NKG7', 'FASLG'],
    'ANTIGEN_PRESENTATION': ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRA', 'HLA-DRB1', 'B2M', 'TAP1', 'TAP2'],
    'CHEMOKINE_SIGNALING': ['CCL2', 'CCL3', 'CCL4', 'CCL5', 'CXCL8', 'CXCL10', 'CCR7', 'CXCR4'],
}

# Combine gene sets
all_gene_sets = {**hallmark_sets, **cell_type_sets}

# Score cells for each gene set
print("Scoring cells for gene set activity...")
pathway_scores = {}
successful_sets = []

for set_name, genes in all_gene_sets.items():
    # Filter to genes present in dataset
    available_genes = [g for g in genes if g in adata.var_names]

    if len(available_genes) >= 2:
        try:
            score_name = f'{set_name}_score'
            sc.tl.score_genes(adata, available_genes, score_name=score_name, use_raw=False)
            pathway_scores[set_name] = adata.obs[score_name].values
            successful_sets.append(set_name)
            print(f"  {set_name}: {len(available_genes)}/{len(genes)} genes")
        except Exception as e:
            print(f"  {set_name}: Failed - {e}")
    else:
        print(f"  {set_name}: Skipped (only {len(available_genes)} genes found)")

print(f"Successfully scored {len(successful_sets)} gene sets")

# Save pathway scores
if pathway_scores:
    scores_df = pd.DataFrame(pathway_scores, index=adata.obs_names)
    scores_df.to_csv('pathway_scores.csv')
else:
    scores_df = pd.DataFrame({'cell': adata.obs_names, 'no_scores': [0] * n_cells})
    scores_df.to_csv('pathway_scores.csv', index=False)

# Perform enrichment analysis on marker genes using gseapy (if available)
print("Performing gene set enrichment analysis...")
enrichment_results = []

try:
    import gseapy as gp

    # Run enrichment for each cluster's top markers
    for cluster in clusters:
        cluster_str = str(cluster)
        cluster_markers = markers_df[markers_df['cluster'] == cluster_str].copy()

        if len(cluster_markers) == 0:
            continue

        # Get top upregulated genes
        top_genes = cluster_markers.nlargest(n_top_genes_param, 'score')['gene'].tolist()

        if len(top_genes) < 3:
            continue

        # Run Enrichr analysis
        try:
            # Try multiple gene set libraries
            libraries = ['GO_Biological_Process_2021', 'KEGG_2021_Human', 'Reactome_2022']

            for library in libraries:
                try:
                    enr = gp.enrichr(
                        gene_list=top_genes,
                        gene_sets=library,
                        organism='Human',
                        outdir=None,
                        no_plot=True,
                        cutoff=0.5
                    )

                    if enr.results is not None and len(enr.results) > 0:
                        results = enr.results.head(10).copy()
                        results['cluster'] = cluster_str
                        results['library'] = library
                        enrichment_results.append(results)
                        print(f"  Cluster {cluster_str} - {library}: {len(results)} terms")
                except Exception as e:
                    print(f"  Cluster {cluster_str} - {library}: Failed ({str(e)[:50]})")

        except Exception as e:
            print(f"  Cluster {cluster_str}: Enrichr failed - {e}")

except ImportError:
    print("WARNING: gseapy not installed. Skipping Enrichr analysis.")
    print("Install with: pip install gseapy")

    # Use simple overlap-based enrichment as fallback
    print("Using built-in pathway overlap analysis instead...")

    for cluster in clusters:
        cluster_str = str(cluster)
        cluster_markers = markers_df[markers_df['cluster'] == cluster_str].copy()

        if len(cluster_markers) == 0:
            continue

        top_genes = set(cluster_markers.nlargest(n_top_genes_param, 'score')['gene'].tolist())

        # Check overlap with our gene sets
        for set_name, genes in all_gene_sets.items():
            overlap = top_genes.intersection(set(genes))
            if len(overlap) >= 2:
                # Simple hypergeometric-like score (overlap / min of both sets)
                overlap_score = len(overlap) / min(len(top_genes), len(genes))
                enrichment_results.append({
                    'cluster': cluster_str,
                    'Term': set_name,
                    'Overlap': f"{len(overlap)}/{len(genes)}",
                    'P-value': max(0.001, 1 - overlap_score),  # Pseudo p-value
                    'Adjusted P-value': max(0.001, 1 - overlap_score),
                    'Genes': ';'.join(overlap),
                    'library': 'Built-in'
                })

except Exception as e:
    print(f"WARNING: Enrichment analysis failed: {e}")

# Save enrichment results
if enrichment_results:
    if isinstance(enrichment_results[0], pd.DataFrame):
        enrichment_df = pd.concat(enrichment_results, ignore_index=True)
    else:
        enrichment_df = pd.DataFrame(enrichment_results)
    enrichment_df.to_csv('enrichment_results.csv', index=False)
    print(f"Saved enrichment results: {len(enrichment_df)} terms")
else:
    enrichment_df = pd.DataFrame({'cluster': [], 'Term': [], 'P-value': []})
    enrichment_df.to_csv('enrichment_results.csv', index=False)
    print("No enrichment results generated")

# Create plots
print("Generating GSEA plots...")
fig = plt.figure(figsize=(16, 12))

# Plot 1: Heatmap of pathway scores by cluster
ax1 = plt.subplot(2, 2, 1)
if len(successful_sets) > 0:
    # Calculate mean pathway score per cluster
    cluster_means = []
    for cluster in clusters:
        cluster_mask = adata.obs[cluster_key] == cluster
        means = {}
        for set_name in successful_sets[:15]:  # Top 15 pathways
            score_col = f'{set_name}_score'
            if score_col in adata.obs.columns:
                means[set_name] = adata.obs[cluster_mask][score_col].mean()
        cluster_means.append(means)

    if cluster_means and cluster_means[0]:
        heatmap_df = pd.DataFrame(cluster_means, index=[str(c) for c in clusters])
        # Z-score normalize for visualization
        heatmap_norm = (heatmap_df - heatmap_df.mean()) / (heatmap_df.std() + 1e-10)

        im = ax1.imshow(heatmap_norm.T.values, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
        ax1.set_xticks(range(len(clusters)))
        ax1.set_xticklabels([str(c) for c in clusters], rotation=45, ha='right')
        ax1.set_yticks(range(len(heatmap_df.columns)))
        ax1.set_yticklabels([s.replace('HALLMARK_', '').replace('_', ' ')[:25] for s in heatmap_df.columns], fontsize=7)
        ax1.set_xlabel('Cluster')
        ax1.set_ylabel('Pathway')
        ax1.set_title('Pathway Activity by Cluster (Z-score)')
        plt.colorbar(im, ax=ax1, label='Z-score')
else:
    ax1.text(0.5, 0.5, 'No pathway scores available', ha='center', va='center')
    ax1.axis('off')

# Plot 2: Top enriched terms (if available)
ax2 = plt.subplot(2, 2, 2)
if len(enrichment_df) > 0 and 'P-value' in enrichment_df.columns:
    # Get top terms by significance
    top_terms = enrichment_df.nsmallest(10, 'P-value')[['Term', 'P-value']].drop_duplicates('Term')

    if len(top_terms) > 0:
        y_pos = range(len(top_terms))
        bars = ax2.barh(y_pos, -np.log10(top_terms['P-value'].values), color='coral', alpha=0.7)
        ax2.set_yticks(y_pos)
        ax2.set_yticklabels([t[:40] for t in top_terms['Term'].values], fontsize=8)
        ax2.set_xlabel('-log10(P-value)')
        ax2.set_title('Top Enriched Terms')
        ax2.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='p=0.05')
        ax2.legend(loc='lower right')
else:
    ax2.text(0.5, 0.5, 'No enrichment results available', ha='center', va='center')
    ax2.axis('off')

# Plot 3: UMAP colored by top pathway score
ax3 = plt.subplot(2, 2, 3)
if len(successful_sets) > 0 and 'X_umap' in adata.obsm:
    top_pathway = successful_sets[0]
    score_col = f'{top_pathway}_score'
    sc.pl.umap(adata, color=score_col, ax=ax3, show=False,
               title=f'{top_pathway.replace("HALLMARK_", "").replace("_", " ")} Score')
elif len(successful_sets) > 0:
    top_pathway = successful_sets[0]
    score_col = f'{top_pathway}_score'
    sc.pl.pca(adata, color=score_col, ax=ax3, show=False,
              title=f'{top_pathway.replace("HALLMARK_", "").replace("_", " ")} Score')
else:
    ax3.text(0.5, 0.5, 'No pathway scores to display', ha='center', va='center')
    ax3.axis('off')

# Plot 4: Distribution of pathway scores
ax4 = plt.subplot(2, 2, 4)
if len(successful_sets) >= 3:
    # Show distribution of top 3 pathways
    for i, set_name in enumerate(successful_sets[:3]):
        score_col = f'{set_name}_score'
        ax4.hist(adata.obs[score_col], bins=30, alpha=0.5,
                label=set_name.replace('HALLMARK_', '').replace('_', ' ')[:20])
    ax4.set_xlabel('Pathway Score')
    ax4.set_ylabel('Number of Cells')
    ax4.set_title('Distribution of Pathway Scores')
    ax4.legend(loc='upper right', fontsize=8)
else:
    ax4.text(0.5, 0.5, 'Insufficient pathways scored', ha='center', va='center')
    ax4.axis('off')

plt.suptitle('Gene Set Enrichment Analysis Results', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('gsea_plots.pdf', dpi=150, bbox_inches='tight')
plt.close()

# Save data
print("Saving data with pathway scores...")
adata.write('gsea_results.h5ad')

# Write summary
summary = f"""Gene Set Enrichment Analysis Summary
=======================================
Input cells: {n_cells}
Input genes: {n_genes}
Cluster key: {cluster_key}
Number of clusters: {n_clusters}

Pathway Scoring:
  Gene sets tested: {len(all_gene_sets)}
  Successfully scored: {len(successful_sets)}
  Top genes per cluster: {n_top_genes_param}

"""

if len(successful_sets) > 0:
    summary += "Scored Pathways:\\n"
    for set_name in successful_sets[:10]:
        summary += f"  - {set_name.replace('HALLMARK_', '').replace('_', ' ')}\\n"
    if len(successful_sets) > 10:
        summary += f"  ... and {len(successful_sets) - 10} more\\n"

summary += f"""
Enrichment Analysis:
  Total enriched terms: {len(enrichment_df)}
"""

if len(enrichment_df) > 0 and 'P-value' in enrichment_df.columns:
    top_5 = enrichment_df.nsmallest(5, 'P-value')
    summary += "  Top 5 enriched terms:\\n"
    for _, row in top_5.iterrows():
        summary += f"    - {row['Term'][:50]} (p={row['P-value']:.2e})\\n"

summary += f"""
Output Files:
  - gsea_results.h5ad: AnnData with pathway scores
  - enrichment_results.csv: GO/pathway enrichment results
  - pathway_scores.csv: Cell-level pathway activity scores
  - gsea_plots.pdf: Visualization plots

Pathway scores stored in adata.obs:
  - '*_score': Activity score for each pathway
"""

with open('gsea_summary.txt', 'w') as f:
    f.write(summary)

print(summary)
print("Gene set enrichment analysis complete!")
    '''
}
