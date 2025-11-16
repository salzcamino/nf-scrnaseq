process TRAJECTORY_ANALYSIS {
    tag "trajectory"
    label 'process_medium'
    publishDir "${params.outdir}/trajectory", mode: params.publish_dir_mode

    input:
    path adata
    val root_cell
    val n_dcs
    val cluster_key

    output:
    path "trajectory_results.h5ad", emit: adata
    path "pseudotime.csv", emit: pseudotime
    path "trajectory_plots.pdf", emit: plots
    path "trajectory_summary.txt", emit: summary

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
    root_cell_param = '!{root_cell}'
    n_dcs_param = int('!{n_dcs}')
    cluster_key_param = '!{cluster_key}'

    # Load data
    adata = sc.read_h5ad('!{adata}')

    summary_lines = []
    summary_lines.append("Trajectory Analysis Summary")
    summary_lines.append("=" * 50)
    summary_lines.append(f"Input cells: {adata.n_obs}")
    summary_lines.append(f"Input genes: {adata.n_vars}")
    summary_lines.append(f"Number of diffusion components: {n_dcs_param}")
    summary_lines.append("")

    # Determine cluster key
    if cluster_key_param == 'auto':
        for key in ['leiden', 'louvain', 'seurat_clusters']:
            if key in adata.obs.columns:
                cluster_key_param = key
                break
        if cluster_key_param == 'auto':
            summary_lines.append("WARNING: No clustering found. Using all cells as one group.")
            cluster_key_param = None
    elif cluster_key_param not in adata.obs.columns:
        summary_lines.append(f"WARNING: Cluster key '{cluster_key_param}' not found")
        cluster_key_param = None

    if cluster_key_param:
        summary_lines.append(f"Cluster key: {cluster_key_param}")
        n_clusters = adata.obs[cluster_key_param].nunique()
        summary_lines.append(f"Number of clusters: {n_clusters}")
    summary_lines.append("")

    # Ensure neighbors are computed
    if 'neighbors' not in adata.uns:
        summary_lines.append("Computing neighbors graph...")
        sc.pp.neighbors(adata)

    # ========================================
    # PAGA (Partition-based Graph Abstraction)
    # ========================================
    summary_lines.append("PAGA Analysis")
    summary_lines.append("-" * 50)

    paga_success = False
    if cluster_key_param:
        try:
            sc.tl.paga(adata, groups=cluster_key_param)
            paga_success = True
            summary_lines.append("PAGA graph computed successfully")

            # Get strong connections
            paga_adj = adata.uns['paga']['connectivities'].toarray()
            n_strong_connections = np.sum(paga_adj > 0.3) // 2  # Divide by 2 for undirected
            summary_lines.append(f"Strong connections (>0.3): {n_strong_connections}")

        except Exception as e:
            summary_lines.append(f"PAGA computation failed: {str(e)}")
    else:
        summary_lines.append("PAGA skipped: requires clustering")

    summary_lines.append("")

    # ========================================
    # Diffusion Pseudotime
    # ========================================
    summary_lines.append("Diffusion Pseudotime Analysis")
    summary_lines.append("-" * 50)

    dpt_success = False
    try:
        # Compute diffusion map
        summary_lines.append(f"Computing diffusion map with {n_dcs_param} components...")
        sc.tl.diffmap(adata, n_comps=n_dcs_param)
        summary_lines.append("Diffusion map computed successfully")

        # Determine root cell
        if root_cell_param == 'auto':
            # Auto-select root based on various criteria
            summary_lines.append("Auto-selecting root cell...")

            # Strategy 1: Use cell with lowest first diffusion component
            # (often corresponds to least differentiated state)
            dc1_values = adata.obsm['X_diffmap'][:, 0]
            root_idx = np.argmin(dc1_values)
            root_cell_id = adata.obs_names[root_idx]

            # Alternative: if we have a stem/progenitor marker
            stem_markers = ['CD34', 'KIT', 'PROM1', 'THY1']
            for marker in stem_markers:
                if marker in adata.var_names:
                    # Use cell with highest expression of stem marker
                    if hasattr(adata.X, 'toarray'):
                        marker_expr = adata[:, marker].X.toarray().flatten()
                    else:
                        marker_expr = adata[:, marker].X.flatten()
                    root_idx = np.argmax(marker_expr)
                    root_cell_id = adata.obs_names[root_idx]
                    summary_lines.append(f"Root selected based on {marker} expression")
                    break

            adata.uns['iroot'] = root_idx
            summary_lines.append(f"Root cell index: {root_idx}")
            summary_lines.append(f"Root cell ID: {root_cell_id}")

        elif root_cell_param.isdigit():
            # Root specified as cluster number
            if cluster_key_param:
                cluster_id = root_cell_param
                cluster_cells = adata.obs[cluster_key_param] == cluster_id
                if cluster_cells.sum() > 0:
                    # Use cell in specified cluster with lowest DC1
                    cluster_indices = np.where(cluster_cells)[0]
                    dc1_in_cluster = adata.obsm['X_diffmap'][cluster_indices, 0]
                    local_idx = np.argmin(dc1_in_cluster)
                    root_idx = cluster_indices[local_idx]
                    adata.uns['iroot'] = root_idx
                    summary_lines.append(f"Root selected from cluster {cluster_id}")
                    summary_lines.append(f"Root cell index: {root_idx}")
                else:
                    summary_lines.append(f"WARNING: Cluster {cluster_id} not found")
                    root_idx = 0
                    adata.uns['iroot'] = root_idx
            else:
                root_idx = int(root_cell_param)
                adata.uns['iroot'] = root_idx
                summary_lines.append(f"Root cell index: {root_idx}")
        else:
            # Root specified as cell ID
            if root_cell_param in adata.obs_names:
                root_idx = adata.obs_names.get_loc(root_cell_param)
                adata.uns['iroot'] = root_idx
                summary_lines.append(f"Root cell: {root_cell_param}")
            else:
                summary_lines.append(f"WARNING: Root cell '{root_cell_param}' not found")
                root_idx = 0
                adata.uns['iroot'] = root_idx

        # Compute diffusion pseudotime
        summary_lines.append("\\nComputing diffusion pseudotime...")
        sc.tl.dpt(adata)
        dpt_success = True

        # DPT statistics
        dpt_values = adata.obs['dpt_pseudotime']
        summary_lines.append(f"Pseudotime range: {dpt_values.min():.4f} - {dpt_values.max():.4f}")
        summary_lines.append(f"Mean pseudotime: {dpt_values.mean():.4f}")
        summary_lines.append(f"Median pseudotime: {dpt_values.median():.4f}")

        # Cells at different pseudotime stages
        early_cells = (dpt_values < 0.33).sum()
        mid_cells = ((dpt_values >= 0.33) & (dpt_values < 0.66)).sum()
        late_cells = (dpt_values >= 0.66).sum()

        summary_lines.append(f"\\nCells by pseudotime stage:")
        summary_lines.append(f"  Early (0-0.33): {early_cells} ({100*early_cells/adata.n_obs:.1f}%)")
        summary_lines.append(f"  Middle (0.33-0.66): {mid_cells} ({100*mid_cells/adata.n_obs:.1f}%)")
        summary_lines.append(f"  Late (0.66-1.0): {late_cells} ({100*late_cells/adata.n_obs:.1f}%)")

        # Pseudotime per cluster
        if cluster_key_param:
            summary_lines.append(f"\\nMean pseudotime per {cluster_key_param} cluster:")
            pt_by_cluster = adata.obs.groupby(cluster_key_param)['dpt_pseudotime'].mean().sort_values()
            for cluster, mean_pt in pt_by_cluster.items():
                summary_lines.append(f"  {cluster}: {mean_pt:.4f}")

    except Exception as e:
        summary_lines.append(f"Diffusion pseudotime failed: {str(e)}")
        import traceback
        summary_lines.append(traceback.format_exc())

    summary_lines.append("")

    # ========================================
    # Force-directed Layout (optional)
    # ========================================
    if paga_success:
        summary_lines.append("Force-directed Layout")
        summary_lines.append("-" * 50)
        try:
            # Initialize with PAGA
            sc.tl.draw_graph(adata, init_pos='paga')
            summary_lines.append("Force-directed layout computed (PAGA-initialized)")
        except Exception as e:
            summary_lines.append(f"Force-directed layout failed: {str(e)}")
            try:
                # Try without PAGA initialization
                sc.tl.draw_graph(adata)
                summary_lines.append("Force-directed layout computed (random initialization)")
            except:
                summary_lines.append("Force-directed layout not available")
        summary_lines.append("")

    # Save pseudotime data
    pt_data = {'cell_id': adata.obs_names}

    if dpt_success:
        pt_data['dpt_pseudotime'] = adata.obs['dpt_pseudotime']
        if 'dpt_groups' in adata.obs.columns:
            pt_data['dpt_groups'] = adata.obs['dpt_groups']

    if 'X_diffmap' in adata.obsm:
        for i in range(min(3, adata.obsm['X_diffmap'].shape[1])):
            pt_data[f'DC{i+1}'] = adata.obsm['X_diffmap'][:, i]

    pt_df = pd.DataFrame(pt_data)
    pt_df.to_csv('pseudotime.csv', index=False)

    # Store trajectory info
    adata.uns['trajectory_analysis'] = {
        'paga_computed': paga_success,
        'dpt_computed': dpt_success,
        'n_diffusion_components': n_dcs_param,
        'cluster_key': cluster_key_param
    }

    # ========================================
    # VISUALIZATION
    # ========================================
    with PdfPages('trajectory_plots.pdf') as pdf:
        # Plot 1: PAGA graph
        if paga_success and cluster_key_param:
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))

            # PAGA connectivity
            ax = axes[0]
            sc.pl.paga(adata, ax=ax, show=False, frameon=False,
                      node_size_scale=1.5, edge_width_scale=1.0)
            ax.set_title(f'PAGA Graph ({cluster_key_param})')

            # PAGA on UMAP
            ax = axes[1]
            if 'X_umap' in adata.obsm:
                sc.pl.paga_compare(adata, ax=ax, show=False, frameon=False,
                                   legend_loc='on data', legend_fontsize=8)
                ax.set_title('PAGA Overlaid on UMAP')
            else:
                ax.text(0.5, 0.5, 'UMAP not available',
                       ha='center', va='center')
                ax.axis('off')

            plt.suptitle('Partition-based Graph Abstraction', fontsize=14)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()

        # Plot 2: Diffusion map
        if 'X_diffmap' in adata.obsm:
            fig, axes = plt.subplots(1, 3, figsize=(18, 5))

            # DC1 vs DC2
            ax = axes[0]
            if dpt_success:
                scatter = ax.scatter(adata.obsm['X_diffmap'][:, 0],
                                    adata.obsm['X_diffmap'][:, 1],
                                    c=adata.obs['dpt_pseudotime'],
                                    cmap='viridis', s=10, alpha=0.7)
                plt.colorbar(scatter, ax=ax, label='Pseudotime')
            else:
                ax.scatter(adata.obsm['X_diffmap'][:, 0],
                          adata.obsm['X_diffmap'][:, 1],
                          s=10, alpha=0.7)
            ax.set_xlabel('DC1')
            ax.set_ylabel('DC2')
            ax.set_title('Diffusion Components 1 & 2')

            # DC1 vs DC3
            ax = axes[1]
            if adata.obsm['X_diffmap'].shape[1] >= 3:
                if dpt_success:
                    scatter = ax.scatter(adata.obsm['X_diffmap'][:, 0],
                                        adata.obsm['X_diffmap'][:, 2],
                                        c=adata.obs['dpt_pseudotime'],
                                        cmap='viridis', s=10, alpha=0.7)
                    plt.colorbar(scatter, ax=ax, label='Pseudotime')
                else:
                    ax.scatter(adata.obsm['X_diffmap'][:, 0],
                              adata.obsm['X_diffmap'][:, 2],
                              s=10, alpha=0.7)
                ax.set_xlabel('DC1')
                ax.set_ylabel('DC3')
                ax.set_title('Diffusion Components 1 & 3')
            else:
                ax.text(0.5, 0.5, 'Only 2 DCs available',
                       ha='center', va='center')
                ax.axis('off')

            # DC2 vs DC3
            ax = axes[2]
            if adata.obsm['X_diffmap'].shape[1] >= 3:
                if cluster_key_param:
                    for cluster in adata.obs[cluster_key_param].unique():
                        mask = adata.obs[cluster_key_param] == cluster
                        ax.scatter(adata.obsm['X_diffmap'][mask, 1],
                                  adata.obsm['X_diffmap'][mask, 2],
                                  s=10, alpha=0.7, label=cluster)
                    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left',
                             title=cluster_key_param, fontsize=8)
                else:
                    ax.scatter(adata.obsm['X_diffmap'][:, 1],
                              adata.obsm['X_diffmap'][:, 2],
                              s=10, alpha=0.7)
                ax.set_xlabel('DC2')
                ax.set_ylabel('DC3')
                ax.set_title('Diffusion Components 2 & 3')
            else:
                ax.text(0.5, 0.5, 'Only 2 DCs available',
                       ha='center', va='center')
                ax.axis('off')

            plt.suptitle('Diffusion Map', fontsize=14)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()

        # Plot 3: Pseudotime on UMAP
        if dpt_success and 'X_umap' in adata.obsm:
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))

            # Pseudotime gradient
            ax = axes[0]
            scatter = ax.scatter(adata.obsm['X_umap'][:, 0],
                                adata.obsm['X_umap'][:, 1],
                                c=adata.obs['dpt_pseudotime'],
                                cmap='viridis', s=10, alpha=0.7)
            plt.colorbar(scatter, ax=ax, label='Pseudotime')
            ax.set_xlabel('UMAP1')
            ax.set_ylabel('UMAP2')
            ax.set_title('Diffusion Pseudotime on UMAP')

            # Mark root cell
            if 'iroot' in adata.uns:
                root_idx = adata.uns['iroot']
                ax.scatter(adata.obsm['X_umap'][root_idx, 0],
                          adata.obsm['X_umap'][root_idx, 1],
                          c='red', s=200, marker='*', edgecolors='black',
                          linewidth=1.5, zorder=5, label='Root')
                ax.legend()

            # Pseudotime distribution
            ax = axes[1]
            ax.hist(adata.obs['dpt_pseudotime'], bins=50, edgecolor='black')
            ax.set_xlabel('Pseudotime')
            ax.set_ylabel('Number of Cells')
            ax.set_title('Pseudotime Distribution')
            ax.axvline(x=adata.obs['dpt_pseudotime'].mean(), color='red',
                      linestyle='--', label=f"Mean: {adata.obs['dpt_pseudotime'].mean():.3f}")
            ax.legend()

            plt.suptitle('Diffusion Pseudotime Analysis', fontsize=14)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()

        # Plot 4: Pseudotime per cluster
        if dpt_success and cluster_key_param:
            fig, ax = plt.subplots(figsize=(10, 6))

            # Violin plot of pseudotime per cluster
            clusters = sorted(adata.obs[cluster_key_param].unique())
            pt_data_plot = [adata.obs.loc[adata.obs[cluster_key_param] == c, 'dpt_pseudotime']
                           for c in clusters]

            parts = ax.violinplot(pt_data_plot, positions=range(len(clusters)),
                                 showmeans=True, showmedians=True)

            # Color violins
            for pc in parts['bodies']:
                pc.set_alpha(0.7)

            ax.set_xticks(range(len(clusters)))
            ax.set_xticklabels(clusters)
            ax.set_xlabel(cluster_key_param.capitalize())
            ax.set_ylabel('Pseudotime')
            ax.set_title(f'Pseudotime Distribution by {cluster_key_param.capitalize()} Cluster')

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()

        # Plot 5: Force-directed layout (if available)
        if 'X_draw_graph_fr' in adata.obsm:
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))

            # By cluster
            ax = axes[0]
            if cluster_key_param:
                for cluster in adata.obs[cluster_key_param].unique():
                    mask = adata.obs[cluster_key_param] == cluster
                    ax.scatter(adata.obsm['X_draw_graph_fr'][mask, 0],
                              adata.obsm['X_draw_graph_fr'][mask, 1],
                              s=10, alpha=0.7, label=cluster)
                ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left',
                         title=cluster_key_param, fontsize=8)
            else:
                ax.scatter(adata.obsm['X_draw_graph_fr'][:, 0],
                          adata.obsm['X_draw_graph_fr'][:, 1],
                          s=10, alpha=0.7)
            ax.set_xlabel('FA1')
            ax.set_ylabel('FA2')
            ax.set_title('Force-Directed Layout (Clusters)')

            # By pseudotime
            ax = axes[1]
            if dpt_success:
                scatter = ax.scatter(adata.obsm['X_draw_graph_fr'][:, 0],
                                    adata.obsm['X_draw_graph_fr'][:, 1],
                                    c=adata.obs['dpt_pseudotime'],
                                    cmap='viridis', s=10, alpha=0.7)
                plt.colorbar(scatter, ax=ax, label='Pseudotime')
            else:
                ax.scatter(adata.obsm['X_draw_graph_fr'][:, 0],
                          adata.obsm['X_draw_graph_fr'][:, 1],
                          s=10, alpha=0.7)
            ax.set_xlabel('FA1')
            ax.set_ylabel('FA2')
            ax.set_title('Force-Directed Layout (Pseudotime)')

            plt.suptitle('PAGA-Initialized Force-Directed Layout', fontsize=14)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()

    # Save results
    adata.write('trajectory_results.h5ad')

    summary_lines.append("Output Files")
    summary_lines.append("-" * 50)
    summary_lines.append("- trajectory_results.h5ad: AnnData with trajectory analysis results")
    summary_lines.append("- pseudotime.csv: Pseudotime values and diffusion components")
    summary_lines.append("- trajectory_plots.pdf: Trajectory visualization plots")
    summary_lines.append("- trajectory_summary.txt: This summary file")

    with open('trajectory_summary.txt', 'w') as f:
        f.write("\\n".join(summary_lines))

    print(f"Trajectory analysis complete. DPT: {dpt_success}, PAGA: {paga_success}")
    '''
}
