process CELL_COMMUNICATION {
    tag "cell_communication"
    label 'process_medium'
    publishDir "${params.outdir}/cell_communication", mode: params.publish_dir_mode

    input:
    path adata
    val cell_type_key
    val min_cells_per_type
    val n_permutations

    output:
    path "communication_results.h5ad", emit: adata
    path "ligand_receptor_pairs.csv", emit: interactions
    path "communication_matrix.csv", emit: matrix
    path "communication_plots.pdf", emit: plots
    path "communication_summary.txt", emit: summary

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
    cell_type_key_param = '!{cell_type_key}'
    min_cells_param = int('!{min_cells_per_type}')
    n_perms_param = int('!{n_permutations}')

    # Load data
    adata = sc.read_h5ad('!{adata}')

    summary_lines = []
    summary_lines.append("Cell-Cell Communication Analysis Summary")
    summary_lines.append("=" * 50)
    summary_lines.append(f"Input cells: {adata.n_obs}")
    summary_lines.append(f"Input genes: {adata.n_vars}")
    summary_lines.append(f"Minimum cells per type: {min_cells_param}")
    summary_lines.append(f"Number of permutations: {n_perms_param}")
    summary_lines.append("")

    # Determine cell type key
    if cell_type_key_param == 'auto':
        # Try common cell type annotation keys
        for key in ['cell_type', 'predicted_labels', 'majority_voting', 'leiden', 'louvain']:
            if key in adata.obs.columns:
                cell_type_key_param = key
                break
        if cell_type_key_param == 'auto':
            summary_lines.append("ERROR: No cell type annotation found")
            summary_lines.append("Please run cell type annotation first or specify cell_type_key")

            # Create minimal outputs
            adata.write('communication_results.h5ad')
            pd.DataFrame().to_csv('ligand_receptor_pairs.csv', index=False)
            pd.DataFrame().to_csv('communication_matrix.csv', index=False)

            with PdfPages('communication_plots.pdf') as pdf:
                fig, ax = plt.subplots()
                ax.text(0.5, 0.5, 'No cell type annotations found', ha='center', va='center')
                ax.axis('off')
                pdf.savefig(fig)
                plt.close()

            with open('communication_summary.txt', 'w') as f:
                f.write("\\n".join(summary_lines))

            import sys
            sys.exit(0)

    summary_lines.append(f"Cell type key: {cell_type_key_param}")

    # Get cell types
    cell_types = adata.obs[cell_type_key_param].astype(str)
    unique_types = cell_types.value_counts()

    summary_lines.append(f"Number of cell types: {len(unique_types)}")
    summary_lines.append("\\nCell type distribution:")
    for ct, count in unique_types.items():
        if count >= min_cells_param:
            summary_lines.append(f"  {ct}: {count} cells")
        else:
            summary_lines.append(f"  {ct}: {count} cells (excluded, < {min_cells_param})")

    # Filter to cell types with enough cells
    valid_types = unique_types[unique_types >= min_cells_param].index.tolist()
    summary_lines.append(f"\\nValid cell types for analysis: {len(valid_types)}")

    if len(valid_types) < 2:
        summary_lines.append("ERROR: Need at least 2 cell types with sufficient cells")
        adata.write('communication_results.h5ad')
        pd.DataFrame().to_csv('ligand_receptor_pairs.csv', index=False)
        pd.DataFrame().to_csv('communication_matrix.csv', index=False)

        with PdfPages('communication_plots.pdf') as pdf:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, 'Insufficient cell types for analysis', ha='center', va='center')
            ax.axis('off')
            pdf.savefig(fig)
            plt.close()

        with open('communication_summary.txt', 'w') as f:
            f.write("\\n".join(summary_lines))

        import sys
        sys.exit(0)

    summary_lines.append("")

    # ========================================
    # TRY CELLPHONEDB FIRST
    # ========================================
    cpdb_success = False
    try:
        from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
        import os
        import tempfile

        summary_lines.append("CellPhoneDB Analysis")
        summary_lines.append("-" * 50)
        summary_lines.append("Using CellPhoneDB for L-R interaction analysis")

        # Prepare data for CellPhoneDB
        # CellPhoneDB expects normalized counts and metadata

        # Create temporary directory for CellPhoneDB outputs
        cpdb_dir = tempfile.mkdtemp()

        # Prepare counts matrix (cells x genes)
        if hasattr(adata.X, 'toarray'):
            counts_df = pd.DataFrame(
                adata.X.toarray().T,
                index=adata.var_names,
                columns=adata.obs_names
            )
        else:
            counts_df = pd.DataFrame(
                adata.X.T,
                index=adata.var_names,
                columns=adata.obs_names
            )

        # Prepare metadata
        meta_df = pd.DataFrame({
            'Cell': adata.obs_names,
            'cell_type': adata.obs[cell_type_key_param].astype(str).values
        })
        meta_df = meta_df.set_index('Cell')

        # Filter to valid cell types
        meta_df = meta_df[meta_df['cell_type'].isin(valid_types)]
        counts_df = counts_df[meta_df.index]

        # Save to temporary files
        counts_path = os.path.join(cpdb_dir, 'counts.txt')
        meta_path = os.path.join(cpdb_dir, 'meta.txt')

        counts_df.to_csv(counts_path, sep='\\t')
        meta_df.to_csv(meta_path, sep='\\t')

        summary_lines.append(f"Cells for analysis: {len(meta_df)}")
        summary_lines.append(f"Genes: {len(counts_df)}")

        # Run CellPhoneDB statistical analysis
        summary_lines.append(f"Running statistical analysis with {n_perms_param} permutations...")

        deconvoluted, means, pvalues, significant_means = cpdb_statistical_analysis_method.call(
            cpdb_file_path=None,  # Use default database
            meta_file_path=meta_path,
            counts_file_path=counts_path,
            counts_data='hgnc_symbol',
            output_path=cpdb_dir,
            iterations=n_perms_param,
            threshold=0.1,
            threads=4,
            debug_seed=-1,
            result_precision=3,
            pvalue=0.05,
            subsampling=False,
            subsampling_log=False,
            subsampling_num_pc=100,
            subsampling_num_cells=None
        )

        summary_lines.append("CellPhoneDB analysis completed successfully!")

        # Process CellPhoneDB results
        # means contains mean expression of L-R pairs
        # pvalues contains p-values from permutation test
        # significant_means contains means where p < threshold

        # Convert to our standard format
        interaction_results = []

        # Get interaction columns (cell type pairs)
        interaction_cols = [col for col in means.columns if '|' in col]

        for idx, row in means.iterrows():
            pair_id = row['id_cp_interaction'] if 'id_cp_interaction' in row else str(idx)
            ligand = row.get('gene_a', row.get('partner_a', 'Unknown'))
            receptor = row.get('gene_b', row.get('partner_b', 'Unknown'))

            for col in interaction_cols:
                sender, receiver = col.split('|')
                if sender in valid_types and receiver in valid_types:
                    mean_val = row[col]
                    if pd.notna(mean_val) and mean_val > 0:
                        # Get p-value
                        p_val = pvalues.loc[idx, col] if col in pvalues.columns else np.nan

                        interaction_results.append({
                            'pair_name': f"{ligand}_{receptor}",
                            'ligand': str(ligand),
                            'receptor': str(receptor),
                            'sender': sender,
                            'receiver': receiver,
                            'ligand_expr': float(mean_val),
                            'receptor_expr': float(mean_val),
                            'communication_score': float(mean_val),
                            'p_value': float(p_val) if pd.notna(p_val) else np.nan,
                            'cpdb_id': pair_id
                        })

        interactions_df = pd.DataFrame(interaction_results)

        if len(interactions_df) > 0:
            interactions_df = interactions_df.sort_values('communication_score', ascending=False)
            interactions_df['significant'] = interactions_df['p_value'] < 0.05

            # Build communication matrix
            communication_matrix = pd.DataFrame(
                0.0,
                index=valid_types,
                columns=valid_types
            )

            for _, row in interactions_df.iterrows():
                communication_matrix.loc[row['sender'], row['receiver']] += row['communication_score']

            n_sig = interactions_df['significant'].sum()
            summary_lines.append(f"\\nTotal interactions found: {len(interactions_df)}")
            summary_lines.append(f"Significant interactions (p<0.05): {n_sig}")
            summary_lines.append(f"Unique L-R pairs: {interactions_df['pair_name'].nunique()}")

        cpdb_success = True

        # Cleanup temp directory
        import shutil
        shutil.rmtree(cpdb_dir, ignore_errors=True)

    except ImportError:
        summary_lines.append("CellPhoneDB not installed - using built-in L-R database")
        summary_lines.append("To use CellPhoneDB: pip install cellphonedb")
        summary_lines.append("")
    except Exception as e:
        summary_lines.append(f"CellPhoneDB analysis failed: {str(e)}")
        summary_lines.append("Falling back to built-in L-R database")
        summary_lines.append("")

    # ========================================
    # FALLBACK: BUILT-IN LIGAND-RECEPTOR DATABASE
    # ========================================
    if not cpdb_success:
        summary_lines.append("Ligand-Receptor Database (Built-in)")
        summary_lines.append("-" * 50)

        # Curated human ligand-receptor pairs (subset of CellPhoneDB/CellChat)
        lr_pairs = {
            # Chemokines and receptors
            'CCL2_CCR2': ('CCL2', 'CCR2'),
            'CCL3_CCR1': ('CCL3', 'CCR1'),
            'CCL4_CCR5': ('CCL4', 'CCR5'),
            'CCL5_CCR5': ('CCL5', 'CCR5'),
            'CXCL8_CXCR1': ('CXCL8', 'CXCR1'),
            'CXCL8_CXCR2': ('CXCL8', 'CXCR2'),
            'CXCL10_CXCR3': ('CXCL10', 'CXCR3'),
            'CXCL12_CXCR4': ('CXCL12', 'CXCR4'),

            # Cytokines
            'IL1B_IL1R1': ('IL1B', 'IL1R1'),
            'IL2_IL2RA': ('IL2', 'IL2RA'),
            'IL4_IL4R': ('IL4', 'IL4R'),
            'IL6_IL6R': ('IL6', 'IL6R'),
            'IL10_IL10RA': ('IL10', 'IL10RA'),
            'IFNG_IFNGR1': ('IFNG', 'IFNGR1'),
            'TNF_TNFRSF1A': ('TNF', 'TNFRSF1A'),
            'TGFB1_TGFBR1': ('TGFB1', 'TGFBR1'),

            # Growth factors
            'EGF_EGFR': ('EGF', 'EGFR'),
            'VEGFA_FLT1': ('VEGFA', 'FLT1'),
            'VEGFA_KDR': ('VEGFA', 'KDR'),
            'FGF2_FGFR1': ('FGF2', 'FGFR1'),
            'PDGFB_PDGFRB': ('PDGFB', 'PDGFRB'),
            'HGF_MET': ('HGF', 'MET'),

            # Immune checkpoints
            'CD274_PDCD1': ('CD274', 'PDCD1'),  # PD-L1/PD-1
            'CD80_CD28': ('CD80', 'CD28'),
            'CD80_CTLA4': ('CD80', 'CTLA4'),
            'CD86_CD28': ('CD86', 'CD28'),
            'CD86_CTLA4': ('CD86', 'CTLA4'),
            'TNFSF9_TNFRSF9': ('TNFSF9', 'TNFRSF9'),  # 4-1BBL/4-1BB

            # MHC interactions
            'HLA-A_CD8A': ('HLA-A', 'CD8A'),
            'HLA-B_CD8A': ('HLA-B', 'CD8A'),
            'HLA-C_CD8A': ('HLA-C', 'CD8A'),
            'HLA-DRA_CD4': ('HLA-DRA', 'CD4'),
            'HLA-DRB1_CD4': ('HLA-DRB1', 'CD4'),

            # Adhesion molecules
            'ICAM1_ITGAL': ('ICAM1', 'ITGAL'),  # ICAM1/LFA-1
            'VCAM1_ITGA4': ('VCAM1', 'ITGA4'),
            'SELP_SELPLG': ('SELP', 'SELPLG'),
            'CDH1_CDH1': ('CDH1', 'CDH1'),  # E-cadherin homotypic

            # Notch signaling
            'JAG1_NOTCH1': ('JAG1', 'NOTCH1'),
            'JAG1_NOTCH2': ('JAG1', 'NOTCH2'),
            'DLL1_NOTCH1': ('DLL1', 'NOTCH1'),
            'DLL4_NOTCH1': ('DLL4', 'NOTCH1'),

            # Wnt signaling
            'WNT5A_FZD5': ('WNT5A', 'FZD5'),
            'WNT3A_FZD1': ('WNT3A', 'FZD1'),

            # Other important pairs
            'CD40LG_CD40': ('CD40LG', 'CD40'),
            'FASLG_FAS': ('FASLG', 'FAS'),
            'TNFSF10_TNFRSF10A': ('TNFSF10', 'TNFRSF10A'),  # TRAIL
            'SPP1_CD44': ('SPP1', 'CD44'),
            'SPP1_ITGAV': ('SPP1', 'ITGAV'),
            'LGALS9_HAVCR2': ('LGALS9', 'HAVCR2'),  # Galectin-9/TIM-3
        }

        # Filter to pairs where both genes are in the dataset
        available_pairs = {}
        for pair_name, (ligand, receptor) in lr_pairs.items():
            if ligand in adata.var_names and receptor in adata.var_names:
                available_pairs[pair_name] = (ligand, receptor)

        summary_lines.append(f"Total L-R pairs in database: {len(lr_pairs)}")
        summary_lines.append(f"Pairs with both genes present: {len(available_pairs)}")

        if len(available_pairs) == 0:
            summary_lines.append("\\nWARNING: No ligand-receptor pairs found in dataset")
            summary_lines.append("This may indicate non-human data or limited gene coverage")

            adata.write('communication_results.h5ad')
            pd.DataFrame().to_csv('ligand_receptor_pairs.csv', index=False)
            pd.DataFrame().to_csv('communication_matrix.csv', index=False)

            with PdfPages('communication_plots.pdf') as pdf:
                fig, ax = plt.subplots()
                ax.text(0.5, 0.5, 'No L-R pairs found in dataset', ha='center', va='center')
                ax.axis('off')
                pdf.savefig(fig)
                plt.close()

            with open('communication_summary.txt', 'w') as f:
                f.write("\\n".join(summary_lines))

            import sys
            sys.exit(0)

        summary_lines.append("\\nAvailable pairs:")
        for pair_name in list(available_pairs.keys())[:10]:
            ligand, receptor = available_pairs[pair_name]
            summary_lines.append(f"  {pair_name}: {ligand} -> {receptor}")
        if len(available_pairs) > 10:
            summary_lines.append(f"  ... and {len(available_pairs) - 10} more")

        summary_lines.append("")

        # ========================================
        # CALCULATE COMMUNICATION SCORES
        # ========================================
        summary_lines.append("Communication Score Calculation")
        summary_lines.append("-" * 50)

        # Get expression matrix (use raw if available, otherwise use .X)
        if adata.raw is not None:
            expr_matrix = adata.raw.to_adata().to_df()
        else:
            if hasattr(adata.X, 'toarray'):
                expr_matrix = pd.DataFrame(
                    adata.X.toarray(),
                    index=adata.obs_names,
                    columns=adata.var_names
                )
            else:
                expr_matrix = pd.DataFrame(
                    adata.X,
                    index=adata.obs_names,
                    columns=adata.var_names
                )

        # Calculate mean expression per cell type
        mean_expr = expr_matrix.groupby(cell_types).mean()
        mean_expr = mean_expr.loc[valid_types]

        # Calculate communication scores for each L-R pair between all cell type pairs
        interaction_results = []
        communication_matrix = pd.DataFrame(
            0.0,
            index=valid_types,
            columns=valid_types
        )

        for pair_name, (ligand, receptor) in available_pairs.items():
            ligand_expr = mean_expr[ligand]
            receptor_expr = mean_expr[receptor]

            for sender in valid_types:
                for receiver in valid_types:
                    # Communication score = ligand_expr(sender) * receptor_expr(receiver)
                    score = float(ligand_expr[sender] * receptor_expr[receiver])

                    if score > 0:
                        # Permutation test for significance
                        if n_perms_param > 0:
                            perm_scores = []
                            sender_cells = expr_matrix.loc[cell_types == sender]
                            receiver_cells = expr_matrix.loc[cell_types == receiver]

                            for _ in range(n_perms_param):
                                # Shuffle cell type labels
                                perm_ligand = np.mean(np.random.choice(
                                    expr_matrix[ligand].values,
                                    size=len(sender_cells)
                                ))
                                perm_receptor = np.mean(np.random.choice(
                                    expr_matrix[receptor].values,
                                    size=len(receiver_cells)
                                ))
                                perm_scores.append(perm_ligand * perm_receptor)

                            p_value = (np.sum(np.array(perm_scores) >= score) + 1) / (n_perms_param + 1)
                        else:
                            p_value = np.nan

                        interaction_results.append({
                            'pair_name': pair_name,
                            'ligand': ligand,
                            'receptor': receptor,
                            'sender': sender,
                            'receiver': receiver,
                            'ligand_expr': float(ligand_expr[sender]),
                            'receptor_expr': float(receptor_expr[receiver]),
                            'communication_score': score,
                            'p_value': p_value
                        })

                        # Add to communication matrix (sum of all interactions)
                        communication_matrix.loc[sender, receiver] += score

        # Create results DataFrame
        interactions_df = pd.DataFrame(interaction_results)

        if len(interactions_df) > 0:
            # Sort by score
            interactions_df = interactions_df.sort_values('communication_score', ascending=False)

            # Add significance flag
            if n_perms_param > 0:
                interactions_df['significant'] = interactions_df['p_value'] < 0.05

                n_sig = interactions_df['significant'].sum()
                summary_lines.append(f"Total interactions scored: {len(interactions_df)}")
                summary_lines.append(f"Significant interactions (p<0.05): {n_sig}")
            else:
                interactions_df['significant'] = True
                summary_lines.append(f"Total interactions scored: {len(interactions_df)}")
                summary_lines.append("Significance testing: disabled (n_permutations=0)")

            # Top interactions
            summary_lines.append("\\nTop 10 interactions by score:")
            for idx, row in interactions_df.head(10).iterrows():
                sig_str = "*" if row.get('significant', True) else ""
                summary_lines.append(f"  {row['sender']} -> {row['receiver']}: "
                                    f"{row['ligand']}-{row['receptor']} "
                                    f"(score={row['communication_score']:.4f}){sig_str}")

            # Most active senders
            sender_activity = interactions_df.groupby('sender')['communication_score'].sum().sort_values(ascending=False)
            summary_lines.append("\\nMost active sender cell types:")
            for sender, total_score in sender_activity.head(5).items():
                summary_lines.append(f"  {sender}: {total_score:.4f}")

            # Most active receivers
            receiver_activity = interactions_df.groupby('receiver')['communication_score'].sum().sort_values(ascending=False)
            summary_lines.append("\\nMost active receiver cell types:")
            for receiver, total_score in receiver_activity.head(5).items():
                summary_lines.append(f"  {receiver}: {total_score:.4f}")

        else:
            summary_lines.append("No interactions found")

    summary_lines.append("")

    # Save interaction results
    interactions_df.to_csv('ligand_receptor_pairs.csv', index=False)
    communication_matrix.to_csv('communication_matrix.csv')

    # Store in AnnData
    n_pairs = len(interactions_df['pair_name'].unique()) if len(interactions_df) > 0 else 0
    adata.uns['cell_communication'] = {
        'cell_type_key': cell_type_key_param,
        'n_pairs_tested': n_pairs,
        'n_interactions': len(interactions_df),
        'valid_cell_types': valid_types,
        'method': 'cellphonedb' if cpdb_success else 'built-in'
    }

    # ========================================
    # VISUALIZATION
    # ========================================
    with PdfPages('communication_plots.pdf') as pdf:
        if len(interactions_df) > 0:
            # Plot 1: Communication heatmap
            fig, ax = plt.subplots(figsize=(10, 8))

            # Normalize matrix by row (sender)
            comm_matrix_norm = communication_matrix.div(communication_matrix.sum(axis=1), axis=0)
            comm_matrix_norm = comm_matrix_norm.fillna(0)

            im = ax.imshow(comm_matrix_norm.values, cmap='YlOrRd', aspect='auto')
            plt.colorbar(im, ax=ax, label='Normalized Communication Score')

            ax.set_xticks(range(len(valid_types)))
            ax.set_yticks(range(len(valid_types)))
            ax.set_xticklabels(valid_types, rotation=45, ha='right')
            ax.set_yticklabels(valid_types)
            ax.set_xlabel('Receiver')
            ax.set_ylabel('Sender')
            ax.set_title('Cell-Cell Communication Matrix')

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()

            # Plot 2: Top L-R pairs
            fig, ax = plt.subplots(figsize=(12, 8))

            # Get top pairs overall
            top_pairs = interactions_df.groupby('pair_name')['communication_score'].sum().sort_values(ascending=True).tail(15)

            colors = plt.cm.viridis(np.linspace(0.3, 0.9, len(top_pairs)))
            ax.barh(range(len(top_pairs)), top_pairs.values, color=colors)
            ax.set_yticks(range(len(top_pairs)))
            ax.set_yticklabels(top_pairs.index)
            ax.set_xlabel('Total Communication Score')
            ax.set_title('Top Ligand-Receptor Pairs by Total Activity')

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()

            # Plot 3: Network visualization (simplified)
            fig, ax = plt.subplots(figsize=(12, 12))

            # Position cell types in a circle
            n_types = len(valid_types)
            angles = np.linspace(0, 2 * np.pi, n_types, endpoint=False)
            x_pos = np.cos(angles)
            y_pos = np.sin(angles)

            # Draw edges (interactions)
            # Only show top interactions to avoid clutter
            top_n_edges = min(50, len(interactions_df))
            top_interactions = interactions_df.head(top_n_edges)

            # Create position dict
            pos_dict = {ct: (x_pos[i], y_pos[i]) for i, ct in enumerate(valid_types)}

            # Draw edges with varying thickness
            max_score = top_interactions['communication_score'].max()
            for _, row in top_interactions.iterrows():
                sender_pos = pos_dict[row['sender']]
                receiver_pos = pos_dict[row['receiver']]

                # Curve the edge
                mid_x = (sender_pos[0] + receiver_pos[0]) / 2
                mid_y = (sender_pos[1] + receiver_pos[1]) / 2

                # Add slight curve toward center
                mid_x *= 0.8
                mid_y *= 0.8

                # Line width based on score
                lw = 0.5 + 3 * (row['communication_score'] / max_score)
                alpha = 0.3 + 0.5 * (row['communication_score'] / max_score)

                ax.plot([sender_pos[0], mid_x, receiver_pos[0]],
                       [sender_pos[1], mid_y, receiver_pos[1]],
                       'b-', alpha=alpha, linewidth=lw)

                # Add arrow
                dx = receiver_pos[0] - mid_x
                dy = receiver_pos[1] - mid_y
                ax.arrow(mid_x, mid_y, dx * 0.3, dy * 0.3,
                        head_width=0.05, head_length=0.03, fc='blue', ec='blue', alpha=alpha)

            # Draw nodes
            for i, ct in enumerate(valid_types):
                ax.scatter(x_pos[i], y_pos[i], s=500, c='lightblue',
                          edgecolors='black', linewidth=2, zorder=5)
                ax.text(x_pos[i] * 1.15, y_pos[i] * 1.15, ct,
                       ha='center', va='center', fontsize=10, fontweight='bold')

            ax.set_xlim(-1.5, 1.5)
            ax.set_ylim(-1.5, 1.5)
            ax.set_aspect('equal')
            ax.axis('off')
            ax.set_title(f'Cell-Cell Communication Network\\n(Top {top_n_edges} interactions)', fontsize=14)

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()

            # Plot 4: Sender vs Receiver activity
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))

            # Sender activity
            ax = axes[0]
            sender_activity = interactions_df.groupby('sender')['communication_score'].sum().sort_values(ascending=True)
            ax.barh(range(len(sender_activity)), sender_activity.values, color='coral')
            ax.set_yticks(range(len(sender_activity)))
            ax.set_yticklabels(sender_activity.index)
            ax.set_xlabel('Total Communication Score')
            ax.set_title('Sender Cell Type Activity')

            # Receiver activity
            ax = axes[1]
            receiver_activity = interactions_df.groupby('receiver')['communication_score'].sum().sort_values(ascending=True)
            ax.barh(range(len(receiver_activity)), receiver_activity.values, color='skyblue')
            ax.set_yticks(range(len(receiver_activity)))
            ax.set_yticklabels(receiver_activity.index)
            ax.set_xlabel('Total Communication Score')
            ax.set_title('Receiver Cell Type Activity')

            plt.suptitle('Cell Type Communication Activity', fontsize=14)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()

            # Plot 5: Dotplot of top interactions
            fig, ax = plt.subplots(figsize=(14, 10))

            # Select top interactions for visualization
            top_for_dot = interactions_df.head(30).copy()
            top_for_dot['interaction'] = top_for_dot['sender'] + ' -> ' + top_for_dot['receiver']
            top_for_dot['pair'] = top_for_dot['ligand'] + '-' + top_for_dot['receptor']

            # Pivot for dotplot
            pivot_score = top_for_dot.pivot_table(
                index='interaction',
                columns='pair',
                values='communication_score',
                fill_value=0
            )

            # Only keep pairs that appear
            pivot_score = pivot_score.loc[:, (pivot_score != 0).any(axis=0)]

            if pivot_score.shape[1] > 0:
                # Create dotplot
                for i, interaction in enumerate(pivot_score.index):
                    for j, pair in enumerate(pivot_score.columns):
                        score = pivot_score.loc[interaction, pair]
                        if score > 0:
                            size = 50 + 200 * (score / pivot_score.values.max())
                            ax.scatter(j, i, s=size, c=score, cmap='YlOrRd',
                                      vmin=0, vmax=pivot_score.values.max(), edgecolors='black', linewidth=0.5)

                ax.set_xticks(range(len(pivot_score.columns)))
                ax.set_yticks(range(len(pivot_score.index)))
                ax.set_xticklabels(pivot_score.columns, rotation=90)
                ax.set_yticklabels(pivot_score.index)
                ax.set_xlabel('Ligand-Receptor Pair')
                ax.set_ylabel('Cell Type Interaction')
                ax.set_title('Top Cell-Cell Communication Interactions')

                # Add colorbar
                sm = plt.cm.ScalarMappable(cmap='YlOrRd',
                                          norm=plt.Normalize(vmin=0, vmax=pivot_score.values.max()))
                sm.set_array([])
                plt.colorbar(sm, ax=ax, label='Communication Score')

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()

        else:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, 'No significant interactions found', ha='center', va='center')
            ax.axis('off')
            pdf.savefig(fig)
            plt.close()

    # Save results
    adata.write('communication_results.h5ad')

    summary_lines.append("Output Files")
    summary_lines.append("-" * 50)
    summary_lines.append("- communication_results.h5ad: AnnData with communication metadata")
    summary_lines.append("- ligand_receptor_pairs.csv: All L-R interaction scores")
    summary_lines.append("- communication_matrix.csv: Cell type communication matrix")
    summary_lines.append("- communication_plots.pdf: Communication visualizations")
    summary_lines.append("- communication_summary.txt: This summary file")

    with open('communication_summary.txt', 'w') as f:
        f.write("\\n".join(summary_lines))

    print(f"Cell communication analysis complete. {len(interactions_df)} interactions scored.")
    '''
}
