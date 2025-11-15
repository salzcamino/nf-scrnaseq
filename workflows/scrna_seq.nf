/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { SAMPLESHEET_CHECK     } from '../modules/local/samplesheet_check'
include { FASTQC                } from '../modules/local/fastqc'
include { MULTIQC               } from '../modules/local/multiqc'
include { STAR_GENOMEGENERATE   } from '../modules/local/star_genomegenerate'
include { STAR_SOLO             } from '../modules/local/star_solo'
include { SEURAT_QC             } from '../modules/local/seurat_qc'
include { SEURAT_ANALYSIS       } from '../modules/local/seurat_analysis'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow SCRNA_SEQ {

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // SUBWORKFLOW: Validate and stage input files
    //
    SAMPLESHEET_CHECK(
        file(params.input)
    )
    ch_versions = ch_versions.mix(SAMPLESHEET_CHECK.out.versions)

    // Parse samplesheet
    SAMPLESHEET_CHECK.out.csv
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [:]
            meta.id = row.sample
            meta.protocol = row.protocol ?: params.protocol

            def reads = []
            if (row.fastq_1) reads.add(file(row.fastq_1, checkIfExists: true))
            if (row.fastq_2) reads.add(file(row.fastq_2, checkIfExists: true))

            return [meta, reads]
        }
        .set { ch_reads }

    //
    // MODULE: Run FastQC
    //
    if (!params.skip_fastqc) {
        FASTQC(ch_reads)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    //
    // MODULE: Prepare genome indices
    //
    if (params.star_index) {
        ch_star_index = Channel.value(file(params.star_index))
    } else {
        ch_fasta = params.fasta ? Channel.value(file(params.fasta, checkIfExists: true)) : Channel.empty()
        ch_gtf = params.gtf ? Channel.value(file(params.gtf, checkIfExists: true)) : Channel.empty()

        STAR_GENOMEGENERATE(
            ch_fasta,
            ch_gtf
        )
        ch_star_index = STAR_GENOMEGENERATE.out.index
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    }

    ch_gtf = params.gtf ? Channel.value(file(params.gtf, checkIfExists: true)) : Channel.empty()

    //
    // MODULE: Run STAR Solo for alignment and quantification
    //
    STAR_SOLO(
        ch_reads,
        ch_star_index,
        ch_gtf,
        params.protocol
    )
    ch_multiqc_files = ch_multiqc_files.mix(STAR_SOLO.out.log_final.collect{it[1]})
    ch_versions = ch_versions.mix(STAR_SOLO.out.versions.first())

    //
    // MODULE: Quality control and cell filtering
    //
    SEURAT_QC(
        STAR_SOLO.out.filtered
    )
    ch_multiqc_files = ch_multiqc_files.mix(SEURAT_QC.out.stats.collect{it[1]})
    ch_versions = ch_versions.mix(SEURAT_QC.out.versions.first())

    //
    // MODULE: Downstream analysis
    //
    if (!params.skip_clustering) {
        SEURAT_ANALYSIS(
            SEURAT_QC.out.filtered
        )
        ch_versions = ch_versions.mix(SEURAT_ANALYSIS.out.versions.first())
    }

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        MULTIQC(
            ch_multiqc_files.collect(),
            Channel.empty(),
            Channel.empty(),
            Channel.empty()
        )
        ch_versions = ch_versions.mix(MULTIQC.out.versions)
    }

    emit:
    multiqc_report = params.skip_multiqc ? Channel.empty() : MULTIQC.out.report
    versions       = ch_versions
}

/*
========================================================================================
    THE END
========================================================================================
*/
