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
""".stripIndent()

// Import modules
include { IMPORT_DATA } from './modules/local/import_data.nf'
include { QC_FILTER } from './modules/local/qc_filter.nf'
include { DOUBLET_DECONTAM } from './modules/local/doublet_decontam.nf'

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

    // Doublet detection and decontamination
    if (params.run_doublet_detection) {
        DOUBLET_DECONTAM(
            QC_FILTER.out.adata,
            params.run_scrublet,
            params.run_scdblfinder,
            params.run_decontx,
            params.scrublet_threshold,
            params.expected_doublet_rate
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
