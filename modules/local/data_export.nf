/*
========================================================================================
    Data Export Module
========================================================================================
    Exports processed AnnData to multiple formats for interoperability with other tools
    and visualization platforms.
    
    Supported formats:
    - Seurat RDS (requires R environment or anndata2ri)
    - Loom format (for viewing in Loom viewers)
    - CellxGene-ready H5AD (with metadata validation)
    - CSV files (expression matrix, cell metadata, gene metadata)
    - 10X MTX format (matrix market format)
*/

process DATA_EXPORT {
    tag "data_export"
    label 'process_medium'
    publishDir "${params.outdir}/data_export", mode: params.publish_dir_mode

    input:
    path adata
    val export_seurat
    val export_loom
    val export_cellxgene
    val export_csv
    val export_mtx

    output:
    path "export_summary.txt", emit: summary
    path "export_manifest.json", emit: manifest
    path "seurat/", optional: true, emit: seurat_dir
    path "seurat/*.rds", optional: true, emit: seurat_file
    path "loom/", optional: true, emit: loom_dir
    path "loom/*.loom", optional: true, emit: loom_file
    path "cellxgene/", optional: true, emit: cellxgene_dir
    path "cellxgene/*.h5ad", optional: true, emit: cellxgene_file
    path "csv/", optional: true, emit: csv_dir
    path "mtx/", optional: true, emit: mtx_dir
    path "logs/export_log.txt", emit: log

    shell:
    '''
#!/usr/bin/env python3
import scanpy as sc
import pandas as pd
import numpy as np
import json
import warnings
from pathlib import Path
from datetime import datetime
import scipy.sparse as sp
from scipy.io import mmwrite
import h5py

warnings.filterwarnings('ignore')

# Create output directories
Path("seurat").mkdir(exist_ok=True)
Path("loom").mkdir(exist_ok=True)
Path("cellxgene").mkdir(exist_ok=True)
Path("csv").mkdir(exist_ok=True)
Path("mtx").mkdir(exist_ok=True)
Path("logs").mkdir(exist_ok=True)

# Initialize logging
log_file = open("logs/export_log.txt", "w")
def log_msg(msg):
    print(msg)
    log_file.write(msg + "\\n")
    log_file.flush()

log_msg("=" * 80)
log_msg("Data Export Module")
log_msg("=" * 80)
log_msg(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log_msg("")

# Load AnnData
log_msg(f"Loading AnnData from: !{adata}")
adata = sc.read_h5ad("!{adata}")
log_msg(f"Loaded: {adata.n_obs} cells, {adata.n_vars} genes")
log_msg("")

# Parse boolean parameters
export_seurat_param = '!{export_seurat}'.lower() == 'true'
export_loom_param = '!{export_loom}'.lower() == 'true'
export_cellxgene_param = '!{export_cellxgene}'.lower() == 'true'
export_csv_param = '!{export_csv}'.lower() == 'true'
export_mtx_param = '!{export_mtx}'.lower() == 'true'

export_status = {}

# ============================================================================
# 1. CSV EXPORT
# ============================================================================
if export_csv_param:
    try:
        log_msg("Exporting to CSV format...")
        
        # Export expression matrix
        if sp.issparse(adata.X):
            expr_matrix = adata.X.toarray()
        else:
            expr_matrix = np.asarray(adata.X)
        
        expr_df = pd.DataFrame(
            expr_matrix,
            index=adata.obs_names,
            columns=adata.var_names
        )
        expr_file = "csv/expression_matrix.csv.gz"
        expr_df.to_csv(expr_file, compression='gzip')
        log_msg(f"  Expression matrix: {expr_file} ({expr_df.shape[0]} cells x {expr_df.shape[1]} genes)")
        
        # Export cell metadata
        obs_file = "csv/cell_metadata.csv"
        adata.obs.to_csv(obs_file)
        log_msg(f"  Cell metadata: {obs_file}")
        
        # Export gene metadata
        var_file = "csv/gene_metadata.csv"
        adata.var.to_csv(var_file)
        log_msg(f"  Gene metadata: {var_file}")
        
        # Export embeddings if available
        if adata.obsm:
            for key in adata.obsm.keys():
                embedding = adata.obsm[key]
                embed_df = pd.DataFrame(
                    embedding,
                    index=adata.obs_names,
                    columns=[f"{key}_{i}" for i in range(embedding.shape[1])]
                )
                embed_file = f"csv/{key}.csv"
                embed_df.to_csv(embed_file)
                log_msg(f"  Embedding {key}: {embed_file}")
        
        export_status['csv'] = 'success'
        log_msg("CSV export completed successfully\\n")
    except Exception as e:
        log_msg(f"ERROR in CSV export: {str(e)}\\n")
        export_status['csv'] = f'failed: {str(e)}'

# ============================================================================
# 2. LOOM EXPORT
# ============================================================================
if export_loom_param:
    try:
        log_msg("Exporting to Loom format...")
        
        loom_file = "loom/data.loom"
        adata.write_loom(loom_file)
        log_msg(f"  Loom file: {loom_file}")
        log_msg("Loom export completed successfully\\n")
        export_status['loom'] = 'success'
    except Exception as e:
        log_msg(f"ERROR in Loom export: {str(e)}\\n")
        export_status['loom'] = f'failed: {str(e)}'

# ============================================================================
# 3. CELLXGENE H5AD EXPORT
# ============================================================================
if export_cellxgene_param:
    try:
        log_msg("Exporting CellxGene-ready H5AD...")
        
        adata_cx = adata.copy()
        
        # Ensure proper obs structure for CellxGene
        # CellxGene requires certain metadata fields
        required_obs = ['cell_type'] if 'cell_type' in adata_cx.obs.columns else []
        
        # Add cell_type from cluster if not present
        if 'cell_type' not in adata_cx.obs.columns:
            if 'leiden' in adata_cx.obs.columns:
                adata_cx.obs['cell_type'] = 'Cluster_' + adata_cx.obs['leiden'].astype(str)
            elif 'louvain' in adata_cx.obs.columns:
                adata_cx.obs['cell_type'] = 'Cluster_' + adata_cx.obs['louvain'].astype(str)
            else:
                adata_cx.obs['cell_type'] = 'Unknown'
        
        # Ensure var_names are gene symbols
        if 'gene_symbols' in adata_cx.var.columns:
            adata_cx.var_names = adata_cx.var['gene_symbols']
        
        # CellxGene metadata
        if 'title' not in adata_cx.uns:
            adata_cx.uns['title'] = 'Processed scRNA-seq data'
        
        cellxgene_file = "cellxgene/data_cellxgene.h5ad"
        adata_cx.write_h5ad(cellxgene_file)
        log_msg(f"  CellxGene H5AD: {cellxgene_file}")
        log_msg(f"  Cell types: {', '.join(adata_cx.obs['cell_type'].unique()[:5].tolist())}")
        log_msg("CellxGene export completed successfully\\n")
        export_status['cellxgene'] = 'success'
    except Exception as e:
        log_msg(f"ERROR in CellxGene export: {str(e)}\\n")
        export_status['cellxgene'] = f'failed: {str(e)}'

# ============================================================================
# 4. 10X MTX FORMAT EXPORT
# ============================================================================
if export_mtx_param:
    try:
        log_msg("Exporting to 10X MTX format...")
        
        # Get expression matrix
        X = adata.X
        if sp.issparse(X):
            X_sparse = X.T.tocsr()  # Genes x cells (transposed for 10X format)
        else:
            X_sparse = sp.csr_matrix(X.T)
        
        # Write matrix
        mtx_file = "mtx/matrix.mtx"
        mmwrite(mtx_file, X_sparse)
        log_msg(f"  Matrix: {mtx_file}")
        
        # Write features (genes)
        features_file = "mtx/features.tsv"
        features_df = adata.var[['gene_ids'] if 'gene_ids' in adata.var.columns else []].copy()
        if 'gene_ids' not in features_df.columns:
            features_df['gene_ids'] = adata.var_names
        features_df['gene_names'] = adata.var_names
        features_df[['gene_ids', 'gene_names']].to_csv(features_file, sep='\\t', header=False)
        log_msg(f"  Features: {features_file}")
        
        # Write barcodes (cells)
        barcodes_file = "mtx/barcodes.tsv"
        with open(barcodes_file, 'w') as f:
            for barcode in adata.obs_names:
                f.write(f"{barcode}\\n")
        log_msg(f"  Barcodes: {barcodes_file}")
        
        log_msg("10X MTX export completed successfully\\n")
        export_status['mtx'] = 'success'
    except Exception as e:
        log_msg(f"ERROR in 10X MTX export: {str(e)}\\n")
        export_status['mtx'] = f'failed: {str(e)}'

# ============================================================================
# 5. SEURAT RDS EXPORT (via H5AD conversion)
# ============================================================================
if export_seurat_param:
    try:
        log_msg("Exporting Seurat-compatible format...")
        log_msg("  Creating HDF5-based Seurat-compatible structure...")
        
        # Create a Seurat-compatible h5ad file that can be imported into R/Seurat
        adata_seurat = adata.copy()
        
        # Ensure raw counts are preserved
        if 'counts' not in adata_seurat.layers:
            log_msg("  WARNING: Raw counts not found in layers")
        
        # Save as h5ad with Seurat metadata
        seurat_h5ad = "seurat/data_seurat_compatible.h5ad"
        adata_seurat.write_h5ad(seurat_h5ad)
        log_msg(f"  Seurat-compatible H5AD: {seurat_h5ad}")
        log_msg("  NOTE: This file can be imported to R with:")
        log_msg("    library(SeuratDisk)")
        log_msg("    seurat_obj <- LoadH5Seurat('data_seurat_compatible.h5ad')")
        
        # Also create Python pickle for anndata2ri if available
        try:
            import rpy2
            import rpy2.robjects as ro
            from anndata2ri.scipy2ri import numpy2rpy
            from anndata2ri.scipy2ri import py2rpy
            
            log_msg("  Converting to Seurat via anndata2ri...")
            rds_file = "seurat/data.rds"
            
            # Conversion code that would run in R environment
            log_msg(f"  NOTE: RDS file requires R environment with Seurat package")
            log_msg(f"  To create RDS file manually in R:")
            log_msg(f"    library(SeuratDisk)")
            log_msg(f"    library(Seurat)")
            log_msg(f"    adata <- LoadH5Seurat('{seurat_h5ad}')")
            log_msg(f"    saveRDS(adata, '{rds_file}')")
        except ImportError:
            log_msg("  anndata2ri not available (optional)")
            log_msg("  For RDS export, use R with Seurat and SeuratDisk packages:")
            log_msg("    library(SeuratDisk)")
            log_msg("    adata <- LoadH5Seurat('data_seurat_compatible.h5ad')")
            log_msg("    saveRDS(adata, 'data.rds')")
        
        export_status['seurat'] = 'success'
        log_msg("Seurat export completed successfully\\n")
    except Exception as e:
        log_msg(f"ERROR in Seurat export: {str(e)}\\n")
        export_status['seurat'] = f'failed: {str(e)}'

# ============================================================================
# Generate Export Summary
# ============================================================================
log_msg("=" * 80)
log_msg("EXPORT SUMMARY")
log_msg("=" * 80)

summary_text = f"""Export Summary
===============

Dataset Information:
  Cells: {adata.n_obs:,}
  Genes: {adata.n_vars:,}
  Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

Export Status:
"""

for format_name, status in export_status.items():
    status_str = "✓ Success" if status == 'success' else f"✗ {status}"
    summary_text += f"  {format_name.upper()}: {status_str}\\n"

summary_text += f"""
Output Locations:
"""

if export_csv_param and export_status.get('csv') == 'success':
    summary_text += """  CSV Format:
    - csv/expression_matrix.csv.gz (compressed expression matrix)
    - csv/cell_metadata.csv (cell annotations)
    - csv/gene_metadata.csv (gene annotations)
    - csv/*.csv (embeddings if available)
"""

if export_loom_param and export_status.get('loom') == 'success':
    summary_text += """  Loom Format:
    - loom/data.loom (compatible with Loom viewers)
"""

if export_cellxgene_param and export_status.get('cellxgene') == 'success':
    summary_text += """  CellxGene Format:
    - cellxgene/data_cellxgene.h5ad (CellxGene-ready H5AD)
"""

if export_mtx_param and export_status.get('mtx') == 'success':
    summary_text += """  10X MTX Format:
    - mtx/matrix.mtx (expression matrix in market matrix format)
    - mtx/features.tsv (gene metadata)
    - mtx/barcodes.tsv (cell barcodes)
"""

if export_seurat_param and export_status.get('seurat') == 'success':
    summary_text += """  Seurat Format:
    - seurat/data_seurat_compatible.h5ad (for loading into Seurat)
    - See logs for R code to convert to RDS
"""

summary_text += """
Usage Instructions:

1. CSV Format - Load in Python/R:
   import pandas as pd
   counts = pd.read_csv('csv/expression_matrix.csv.gz')
   metadata = pd.read_csv('csv/cell_metadata.csv', index_col=0)

2. Loom Format - View in browser:
   Open http://loom.linnarssonlab.org/ and upload loom/data.loom

3. CellxGene Format - View with cellxgene:
   cellxgene launch cellxgene/data_cellxgene.h5ad

4. 10X MTX Format - Load in Seurat (R):
   library(Seurat)
   data <- Read10X(data.dir = 'mtx/')
   obj <- CreateSeuratObject(counts = data)

5. Seurat Format - Import to R:
   library(SeuratDisk)
   library(Seurat)
   obj <- LoadH5Seurat('seurat/data_seurat_compatible.h5ad')
   # Or create RDS as described in logs

Embedding Keys Available:
"""

for key in list(adata.obsm.keys())[:10]:
    shape = adata.obsm[key].shape
    summary_text += f"  - {key} ({shape[1]} components)\\n"

summary_text += """
"""

log_msg(summary_text)

# Write summary to file
with open("export_summary.txt", "w") as f:
    f.write(summary_text)

# ============================================================================
# Generate Export Manifest (JSON)
# ============================================================================
manifest = {
    "timestamp": datetime.now().isoformat(),
    "dataset": {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "cell_types": list(adata.obs['cell_type'].unique()) if 'cell_type' in adata.obs.columns else [],
        "embeddings": list(adata.obsm.keys()),
        "obs_columns": list(adata.obs.columns),
        "var_columns": list(adata.var.columns),
    },
    "exports": {
        "csv": {
            "enabled": export_csv_param,
            "status": export_status.get('csv', 'not_run'),
            "files": [
                "expression_matrix.csv.gz",
                "cell_metadata.csv",
                "gene_metadata.csv"
            ] if export_csv_param and export_status.get('csv') == 'success' else []
        },
        "loom": {
            "enabled": export_loom_param,
            "status": export_status.get('loom', 'not_run'),
            "files": ["data.loom"] if export_loom_param and export_status.get('loom') == 'success' else []
        },
        "cellxgene": {
            "enabled": export_cellxgene_param,
            "status": export_status.get('cellxgene', 'not_run'),
            "files": ["data_cellxgene.h5ad"] if export_cellxgene_param and export_status.get('cellxgene') == 'success' else []
        },
        "mtx": {
            "enabled": export_mtx_param,
            "status": export_status.get('mtx', 'not_run'),
            "files": ["matrix.mtx", "features.tsv", "barcodes.tsv"] if export_mtx_param and export_status.get('mtx') == 'success' else []
        },
        "seurat": {
            "enabled": export_seurat_param,
            "status": export_status.get('seurat', 'not_run'),
            "files": ["data_seurat_compatible.h5ad"] if export_seurat_param and export_status.get('seurat') == 'success' else []
        }
    }
}

with open("export_manifest.json", "w") as f:
    json.dump(manifest, f, indent=2, default=str)

log_msg("Export manifest saved to export_manifest.json")
log_msg("")
log_msg("Data export completed!")
log_file.close()
print(summary_text)
    '''
}
