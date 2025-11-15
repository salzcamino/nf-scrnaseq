process STAR_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::star=2.7.10b bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' :
        'biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' }"

    input:
    tuple val(meta), path(reads)
    path  index
    path  gtf

    output:
    tuple val(meta), path('*Aligned.sortedByCoord.out.bam')  , emit: bam
    tuple val(meta), path('*Log.final.out')                  , emit: log_final
    tuple val(meta), path('*Log.out')                        , emit: log_out
    tuple val(meta), path('*Log.progress.out')               , emit: log_progress
    tuple val(meta), path('*SJ.out.tab')                     , emit: sj
    tuple val(meta), path('*ReadsPerGene.out.tab')           , optional:true, emit: reads_per_gene
    path  "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_files = reads instanceof List ? reads.join(',') : reads
    def memory_gb = task.memory ? task.memory.toGiga() - 1 : 30

    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn $read_files \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix ${prefix}. \\
        --limitBAMsortRAM ${memory_gb * 1000000000} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed 's/STAR_//g')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
