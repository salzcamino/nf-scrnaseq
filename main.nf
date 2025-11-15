#!/usr/bin/env nextflow
/*
========================================================================================
    nf-scrnaseq
========================================================================================
    Single-cell RNA-seq Analysis Pipeline
    Github: https://github.com/salzcamino/nf-scrnaseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    PRINT HELP MESSAGE
========================================================================================
*/

def helpMessage() {
    log.info"""
    =========================================
    nf-scrnaseq v${workflow.manifest.version}
    =========================================

    Usage:
      nextflow run main.nf --input samplesheet.csv --outdir results [options]

    Mandatory arguments:
      --input                       Path to input samplesheet (CSV format)
      --outdir                      Output directory for results
      --genome                      Reference genome (e.g., 'GRCh38', 'GRCm39')
      --fasta                       Path to genome FASTA file (if not using --genome)
      --gtf                         Path to genome GTF annotation file

    Optional arguments:
      --protocol                    scRNA-seq protocol (default: '10x_3prime')
                                    Options: 10x_3prime, 10x_5prime, dropseq, smartseq2
      --skip_fastqc                 Skip FastQC step
      --skip_multiqc                Skip MultiQC step
      --aligner                     Alignment tool (default: 'star')
                                    Options: star, salmon

      Cell filtering parameters:
      --min_genes                   Minimum number of genes per cell (default: 200)
      --min_cells                   Minimum number of cells per gene (default: 3)
      --max_mito_pct                Maximum mitochondrial percentage (default: 5)
      --min_counts                  Minimum total counts per cell (default: 1000)

      Analysis parameters:
      --n_hvgs                      Number of highly variable genes (default: 2000)
      --n_pcs                       Number of principal components (default: 50)
      --resolution                  Clustering resolution (default: 0.8)
      --skip_clustering             Skip clustering analysis
      --skip_annotation             Skip cell type annotation

      Resource options:
      --max_cpus                    Maximum number of CPUs (default: 16)
      --max_memory                  Maximum memory (default: 128.GB)
      --max_time                    Maximum time (default: 240.h)

      Other options:
      -profile                      Configuration profile (docker, singularity, conda, test)
      -resume                       Resume previous run
      --help                        Show this help message

    Example:
      nextflow run main.nf \\
        --input samplesheet.csv \\
        --outdir results \\
        --genome GRCh38 \\
        --protocol 10x_3prime \\
        -profile docker
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// Check mandatory parameters
if (!params.input) {
    error "Please provide an input samplesheet with --input"
}

if (!params.outdir) {
    error "Please provide an output directory with --outdir"
}

if (!params.genome && (!params.fasta || !params.gtf)) {
    error "Please provide either --genome or both --fasta and --gtf"
}

/*
========================================================================================
    IMPORT WORKFLOWS
========================================================================================
*/

include { SCRNA_SEQ } from './workflows/scrna_seq'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {

    // Print pipeline information
    log.info """
    =========================================
    nf-scrnaseq v${workflow.manifest.version}
    =========================================
    Input samplesheet : ${params.input}
    Output directory  : ${params.outdir}
    Genome            : ${params.genome ?: 'Custom'}
    Protocol          : ${params.protocol}
    Aligner           : ${params.aligner}
    =========================================
    """.stripIndent()

    // Run main workflow
    SCRNA_SEQ()
}

/*
========================================================================================
    WORKFLOW COMPLETION
========================================================================================
*/

workflow.onComplete {
    log.info """
    =========================================
    Pipeline completed!
    Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Time      : ${workflow.duration}
    Results   : ${params.outdir}
    =========================================
    """.stripIndent()
}

workflow.onError {
    log.error "Pipeline failed with error: ${workflow.errorMessage}"
}
