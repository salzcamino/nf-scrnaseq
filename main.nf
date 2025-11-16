#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================================
    nf-scrnaseq: Single-cell RNA-seq Analysis Pipeline
========================================================================================
    Github: https://github.com/salzcamino/nf-scrnaseq
----------------------------------------------------------------------------------------
*/

// Print help message if requested
def helpMessage() {
    log.info """
    Usage:
      nextflow run main.nf --input <path> [options]

    Required arguments:
      --input              Path to input data (10X directory, h5ad file, or CSV file)

    Input options:
      --input_format       Input format: 'auto', '10x', 'h5ad', 'csv' (default: auto)

    QC filtering options:
      --min_genes          Minimum number of genes per cell (default: 200)
      --min_cells          Minimum number of cells per gene (default: 3)
      --max_genes          Maximum number of genes per cell (default: 2500)
      --max_counts         Maximum counts per cell (default: auto-detect)
      --max_pct_mt         Maximum mitochondrial percentage (default: 5)

    Gene filtering options:
      --exclude_mt         Exclude mitochondrial genes (default: false)
      --exclude_ribo       Exclude ribosomal genes (default: false)

    Doublet detection options:
      --run_doublet_detection  Run doublet detection (default: true)
      --run_scrublet          Run Scrublet doublet detection (default: true)
      --run_scdblfinder       Run scDblFinder (R-based, slower) (default: false)
      --run_decontx           Run DecontX contamination estimation (default: false)
      --scrublet_threshold    Scrublet threshold: 'auto' or numeric (default: auto)
      --expected_doublet_rate Expected doublet rate (default: 0.06)

    Normalization options:
      --target_sum            Target sum for normalization (default: auto/median)
      --log_transform         Apply log1p transformation (default: true)

    Highly variable genes:
      --n_top_genes           Number of top HVGs to select (default: 2000)
      --hvg_flavor            HVG selection method: 'seurat' or 'cell_ranger' (default: seurat)

    Dimensionality reduction:
      --n_pcs                 Number of principal components (default: 50)
      --n_neighbors           Number of neighbors for UMAP (default: 15)
      --run_umap              Compute UMAP embedding (default: true)
      --run_tsne              Compute t-SNE embedding (default: false)

    Clustering:
      --run_leiden            Run Leiden clustering (default: true)
      --run_louvain           Run Louvain clustering (default: false)
      --leiden_resolution     Leiden resolution parameter (default: 1.0)
      --louvain_resolution    Louvain resolution parameter (default: 1.0)
      --run_seurat_clustering Run Seurat SNN clustering (R-based) (default: false)
      --run_celda             Run Celda clustering (R-based) (default: false)
      --seurat_resolution     Seurat resolution parameter (default: 0.8)

    Differential expression:
      --run_diff_expression   Run differential expression analysis (default: true)
      --de_method             Method: 'wilcoxon', 't-test', 'logreg' (default: wilcoxon)
      --de_n_genes            Number of top genes per cluster (default: 25)
      --de_min_fold_change    Minimum fold change (default: 1.5)

    Cell type annotation:
      --run_annotation        Run cell type annotation (default: true)
      --marker_file           Marker gene file: 'default' or path to JSON/CSV (default: default)
      --annotation_method     Scoring method (default: score_genes)

    Output options:
      --outdir             Output directory (default: ./results)
      --publish_dir_mode   Publishing mode: 'copy', 'symlink', 'move' (default: copy)

    Profiles:
      -profile docker      Use Docker containers
      -profile singularity Use Singularity containers
      -profile conda       Use Conda environment
      -profile test        Run with test data

    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// Validate input parameters
if (!params.input) {
    log.error "ERROR: --input parameter is required!"
    helpMessage()
    exit 1
}

// Print parameter summary
log.info """
=======================================================
nf-scrnaseq v${workflow.manifest.version}
=======================================================
Input Data     : ${params.input}
Input Format   : ${params.input_format}
Output Dir     : ${params.outdir}
-------------------------------------------------------
QC Parameters:
  Min genes    : ${params.min_genes}
  Min cells    : ${params.min_cells}
  Max genes    : ${params.max_genes}
  Max counts   : ${params.max_counts ?: 'auto'}
  Max MT %     : ${params.max_pct_mt}
-------------------------------------------------------
Doublet Detection:
  Enabled      : ${params.run_doublet_detection}
  Scrublet     : ${params.run_scrublet}
  scDblFinder  : ${params.run_scdblfinder}
  DecontX      : ${params.run_decontx}
-------------------------------------------------------
Normalization:
  Target sum   : ${params.target_sum}
  Log transform: ${params.log_transform}
-------------------------------------------------------
HVG Selection:
  N top genes  : ${params.n_top_genes}
  Flavor       : ${params.hvg_flavor}
-------------------------------------------------------
Dim Reduction:
  N PCs        : ${params.n_pcs}
  N neighbors  : ${params.n_neighbors}
  Run UMAP     : ${params.run_umap}
  Run t-SNE    : ${params.run_tsne}
-------------------------------------------------------
Clustering:
  Run Leiden   : ${params.run_leiden}
  Run Louvain  : ${params.run_louvain}
  Leiden res   : ${params.leiden_resolution}
  Louvain res  : ${params.louvain_resolution}
  Run Seurat   : ${params.run_seurat_clustering}
  Run Celda    : ${params.run_celda}
  Seurat res   : ${params.seurat_resolution}
-------------------------------------------------------
Diff Expression:
  Enabled      : ${params.run_diff_expression}
  Method       : ${params.de_method}
  N genes      : ${params.de_n_genes}
  Min FC       : ${params.de_min_fold_change}
-------------------------------------------------------
Cell Type Annotation:
  Enabled      : ${params.run_annotation}
  Marker file  : ${params.marker_file}
  Method       : ${params.annotation_method}
-------------------------------------------------------
""".stripIndent()

// Import modules
include { IMPORT_DATA } from './modules/local/import_data.nf'
include { QC_FILTER } from './modules/local/qc_filter.nf'
include { DOUBLET_DECONTAM } from './modules/local/doublet_decontam.nf'
include { NORMALIZE } from './modules/local/normalize.nf'
include { HIGHLY_VARIABLE_GENES } from './modules/local/highly_variable_genes.nf'
include { REDUCE_DIMS } from './modules/local/reduce_dims.nf'
include { CLUSTERING } from './modules/local/clustering.nf'
include { DIFF_EXPRESSION } from './modules/local/diff_expression.nf'
include { CELL_TYPE_ANNOTATION } from './modules/local/cell_type_annotation.nf'

/*
========================================================================================
    Main Workflow
========================================================================================
*/

workflow {
    // Create input channel
    ch_input = Channel.fromPath(params.input, checkIfExists: true)

    // Import data
    IMPORT_DATA(
        ch_input,
        params.input_format
    )

    // QC and filtering
    QC_FILTER(
        IMPORT_DATA.out.adata,
        params.min_genes,
        params.min_cells,
        params.max_genes,
        params.max_counts ?: 'null',
        params.max_pct_mt,
        params.exclude_mt,
        params.exclude_ribo
    )

    // Doublet detection and decontamination (optional)
    if (params.run_doublet_detection) {
        DOUBLET_DECONTAM(
            QC_FILTER.out.adata,
            params.run_scrublet,
            params.run_scdblfinder,
            params.run_decontx,
            params.scrublet_threshold,
            params.expected_doublet_rate
        )
        ch_for_norm = DOUBLET_DECONTAM.out.adata
    } else {
        ch_for_norm = QC_FILTER.out.adata
    }

    // Normalization
    NORMALIZE(
        ch_for_norm,
        params.target_sum,
        params.log_transform
    )

    // Highly variable gene selection
    HIGHLY_VARIABLE_GENES(
        NORMALIZE.out.adata,
        params.n_top_genes,
        params.hvg_min_mean,
        params.hvg_max_mean,
        params.hvg_min_disp,
        params.hvg_flavor
    )

    // Dimensionality reduction
    REDUCE_DIMS(
        HIGHLY_VARIABLE_GENES.out.adata,
        params.n_pcs,
        params.n_neighbors,
        params.run_umap,
        params.run_tsne,
        params.umap_min_dist,
        params.tsne_perplexity
    )

    // Clustering
    CLUSTERING(
        REDUCE_DIMS.out.adata,
        params.run_leiden,
        params.run_louvain,
        params.leiden_resolution,
        params.louvain_resolution,
        params.cluster_key,
        params.run_seurat_clustering,
        params.run_celda,
        params.seurat_resolution,
        params.celda_L,
        params.celda_K
    )

    // Differential expression analysis (optional)
    if (params.run_diff_expression) {
        DIFF_EXPRESSION(
            CLUSTERING.out.adata,
            params.cluster_key,
            params.de_method,
            params.de_n_genes,
            params.de_min_fold_change,
            params.de_min_in_group_fraction,
            params.de_max_out_group_fraction
        )

        // Cell type annotation (optional, depends on DE)
        if (params.run_annotation) {
            CELL_TYPE_ANNOTATION(
                DIFF_EXPRESSION.out.adata,
                params.marker_file,
                params.cluster_key,
                params.annotation_method
            )
        }
    } else if (params.run_annotation) {
        // Run annotation without DE
        CELL_TYPE_ANNOTATION(
            CLUSTERING.out.adata,
            params.marker_file,
            params.cluster_key,
            params.annotation_method
        )
    }
}

/*
========================================================================================
    Workflow Completion
========================================================================================
*/

workflow.onComplete {
    log.info """
    =======================================================
    Pipeline execution summary
    =======================================================
    Completed at : ${workflow.complete}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Work Dir     : ${workflow.workDir}
    Exit status  : ${workflow.exitStatus}
    =======================================================
    """.stripIndent()
}

workflow.onError {
    log.error "Pipeline execution stopped with error: ${workflow.errorMessage}"
}
