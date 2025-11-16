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
      --annotation_method     Method: 'celltypist' or 'marker_scoring' (default: celltypist)
      --celltypist_model      CellTypist model name (default: Immune_All_Low.pkl)
      --marker_file           Marker gene file for marker_scoring (default: default)

    Gene set enrichment:
      --run_gsea              Run gene set enrichment analysis (default: true)
      --gsea_gene_sets        Gene set database: 'default', 'GO', 'KEGG' (default: default)
      --gsea_n_top_genes      Number of top genes for enrichment (default: 100)

    Batch correction:
      --run_batch_correction  Run batch effect detection/correction (default: true)
      --batch_key             Column containing batch info (default: batch)
      --batch_correction_method  Method: 'harmony', 'combat', 'bbknn' (default: harmony)
      --batch_effect_threshold   Threshold for applying correction (default: 0.3)

    Cell cycle scoring:
      --run_cell_cycle        Run cell cycle phase scoring (default: true)
      --regress_cell_cycle    Regress out cell cycle effects (default: false)

    Trajectory analysis:
      --run_trajectory        Run trajectory/pseudotime analysis (default: true)
      --trajectory_root       Root cell: 'auto', cluster ID, or barcode (default: auto)
      --n_diffusion_comps     Number of diffusion components (default: 15)

    Cell-cell communication:
      --run_communication     Run cell-cell communication analysis (default: true)
      --communication_cell_type_key  Cell type column (default: auto)
      --communication_min_cells      Min cells per type (default: 10)
      --communication_n_permutations  Permutations for p-value (default: 100)

    HTML Report:
      --generate_report       Generate HTML report (default: true)
      --report_title          Report title (default: nf-scrnaseq Analysis Report)

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
  Method       : ${params.annotation_method}
  Model        : ${params.celltypist_model}
-------------------------------------------------------
Gene Set Enrichment:
  Enabled      : ${params.run_gsea}
  Gene sets    : ${params.gsea_gene_sets}
  Top genes    : ${params.gsea_n_top_genes}
-------------------------------------------------------
Batch Correction:
  Enabled      : ${params.run_batch_correction}
  Batch key    : ${params.batch_key}
  Method       : ${params.batch_correction_method}
  Threshold    : ${params.batch_effect_threshold}
-------------------------------------------------------
Cell Cycle:
  Enabled      : ${params.run_cell_cycle}
  Regress      : ${params.regress_cell_cycle}
-------------------------------------------------------
Trajectory:
  Enabled      : ${params.run_trajectory}
  Root cell    : ${params.trajectory_root}
  N DCs        : ${params.n_diffusion_comps}
-------------------------------------------------------
Cell Communication:
  Enabled      : ${params.run_communication}
  Cell type key: ${params.communication_cell_type_key}
  Min cells    : ${params.communication_min_cells}
  Permutations : ${params.communication_n_permutations}
-------------------------------------------------------
HTML Report:
  Enabled      : ${params.generate_report}
  Title        : ${params.report_title}
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
include { GENE_SET_ENRICHMENT } from './modules/local/gsea.nf'
include { BATCH_CORRECTION } from './modules/local/batch_correction.nf'
include { CELL_CYCLE_SCORING } from './modules/local/cell_cycle.nf'
include { TRAJECTORY_ANALYSIS } from './modules/local/trajectory.nf'
include { CELL_COMMUNICATION } from './modules/local/cell_communication.nf'
include { HTML_REPORT } from './modules/local/html_report.nf'

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

    // Cell cycle scoring (optional, runs on normalized data)
    if (params.run_cell_cycle) {
        CELL_CYCLE_SCORING(
            NORMALIZE.out.adata,
            params.regress_cell_cycle
        )
        ch_for_hvg = CELL_CYCLE_SCORING.out.adata
    } else {
        ch_for_hvg = NORMALIZE.out.adata
    }

    // Highly variable gene selection
    HIGHLY_VARIABLE_GENES(
        ch_for_hvg,
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

    // Batch correction (optional, runs after dim reduction)
    if (params.run_batch_correction) {
        BATCH_CORRECTION(
            REDUCE_DIMS.out.adata,
            params.batch_key,
            params.batch_correction_method,
            params.batch_effect_threshold
        )
        ch_for_clustering = BATCH_CORRECTION.out.adata
    } else {
        ch_for_clustering = REDUCE_DIMS.out.adata
    }

    // Clustering
    CLUSTERING(
        ch_for_clustering,
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
                params.annotation_method,
                params.celltypist_model,
                params.marker_file,
                params.cluster_key
            )

            // Gene set enrichment analysis (optional, requires DE results)
            if (params.run_gsea) {
                GENE_SET_ENRICHMENT(
                    CELL_TYPE_ANNOTATION.out.adata,
                    DIFF_EXPRESSION.out.markers,
                    params.gsea_gene_sets,
                    params.gsea_n_top_genes,
                    params.cluster_key
                )
            }
        } else if (params.run_gsea) {
            // Run GSEA without annotation
            GENE_SET_ENRICHMENT(
                DIFF_EXPRESSION.out.adata,
                DIFF_EXPRESSION.out.markers,
                params.gsea_gene_sets,
                params.gsea_n_top_genes,
                params.cluster_key
            )
        }
    } else if (params.run_annotation) {
        // Run annotation without DE
        CELL_TYPE_ANNOTATION(
            CLUSTERING.out.adata,
            params.annotation_method,
            params.celltypist_model,
            params.marker_file,
            params.cluster_key
        )
    }

    // Trajectory analysis (optional, runs after clustering)
    if (params.run_trajectory) {
        TRAJECTORY_ANALYSIS(
            CLUSTERING.out.adata,
            params.trajectory_root,
            params.n_diffusion_comps,
            params.cluster_key
        )
    }

    // Cell-cell communication (optional, requires cell type annotations)
    if (params.run_communication) {
        // Use annotated data if available, otherwise clustered data
        if (params.run_annotation) {
            CELL_COMMUNICATION(
                CELL_TYPE_ANNOTATION.out.adata,
                params.communication_cell_type_key,
                params.communication_min_cells,
                params.communication_n_permutations
            )
        } else {
            // Use clusters as cell types
            CELL_COMMUNICATION(
                CLUSTERING.out.adata,
                params.cluster_key,
                params.communication_min_cells,
                params.communication_n_permutations
            )
        }
    }

    // HTML Report Generation (optional, consolidates all results)
    if (params.generate_report) {
        // Determine which AnnData to use for the report (most complete version)
        if (params.run_annotation && params.run_diff_expression) {
            HTML_REPORT(CELL_TYPE_ANNOTATION.out.adata)
        } else if (params.run_diff_expression) {
            HTML_REPORT(DIFF_EXPRESSION.out.adata)
        } else {
            HTML_REPORT(CLUSTERING.out.adata)
        }
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
