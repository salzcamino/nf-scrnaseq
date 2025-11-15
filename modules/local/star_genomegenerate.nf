process STAR_GENOMEGENERATE {
    tag "$fasta"
    label 'process_high'

    conda "bioconda::star=2.7.10b"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.10b--h6b7c446_1' :
        'biocontainers/star:2.7.10b--h6b7c446_1' }"

    input:
    path fasta
    path gtf

    output:
    path "star"         , emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def memory_gb = task.memory ? task.memory.toGiga() : 0
    def avail_mem = memory_gb > 0 ? "--limitGenomeGenerateRAM ${memory_gb * 1024000000}" : ''

    """
    mkdir star

    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star/ \\
        --genomeFastaFiles $fasta \\
        --sjdbGTFfile $gtf \\
        --runThreadN $task.cpus \\
        $avail_mem \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed 's/STAR_//g')
    END_VERSIONS
    """
}
