process STAR_SOLO {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::star=2.7.10b"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.10b--h6b7c446_1' :
        'biocontainers/star:2.7.10b--h6b7c446_1' }"

    input:
    tuple val(meta), path(reads)
    path  index
    path  gtf
    val   protocol

    output:
    tuple val(meta), path('*Solo.out/Gene*/filtered')         , emit: filtered
    tuple val(meta), path('*Solo.out/Gene*/raw')              , emit: raw
    tuple val(meta), path('*Aligned.sortedByCoord.out.bam')   , optional:true, emit: bam
    tuple val(meta), path('*Log.final.out')                   , emit: log_final
    tuple val(meta), path('*Log.out')                         , emit: log_out
    tuple val(meta), path('*SJ.out.tab')                      , emit: sj
    path  "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read1 = reads[0]
    def read2 = reads.size() > 1 ? reads[1] : ''
    def memory_gb = task.memory ? task.memory.toGiga() - 1 : 30

    // Protocol-specific parameters
    def solo_type = 'CB_UMI_Simple'
    def cb_whitelist = ''
    def cb_len = 16
    def umi_len = 12
    def cb_start = 1
    def umi_start = 17

    // Adjust based on protocol
    if (protocol == '10x_3prime') {
        cb_whitelist = '--soloCBwhitelist 737K-august-2016.txt'
        cb_len = 16
        umi_len = 12
    } else if (protocol == '10x_5prime') {
        cb_whitelist = '--soloCBwhitelist 737K-august-2016.txt'
        cb_len = 16
        umi_len = 10
    }

    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn $read2 $read1 \\
        --readFilesCommand zcat \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix ${prefix}. \\
        --outSAMtype BAM SortedByCoordinate \\
        --limitBAMsortRAM ${memory_gb * 1000000000} \\
        --soloType $solo_type \\
        $cb_whitelist \\
        --soloCBstart $cb_start \\
        --soloCBlen $cb_len \\
        --soloUMIstart $umi_start \\
        --soloUMIlen $umi_len \\
        --soloBarcodeReadLength 0 \\
        --soloFeatures Gene GeneFull \\
        --soloUMIdedup 1MM_All \\
        --soloUMIfiltering MultiGeneUMI \\
        --soloCellFilter EmptyDrops_CR \\
        --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed 's/STAR_//g')
    END_VERSIONS
    """
}
