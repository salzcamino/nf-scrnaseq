process SEURAT_QC {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::r-base=4.3.1 bioconda::r-seurat=4.3.0 conda-forge::r-ggplot2=3.4.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-seurat:4.3.0' :
        'satijalab/seurat:4.3.0' }"

    input:
    tuple val(meta), path(matrix_dir)

    output:
    tuple val(meta), path("*_filtered.rds")      , emit: filtered
    tuple val(meta), path("*_qc_*.pdf")          , emit: qc_plots
    tuple val(meta), path("*_filtering_stats.csv"), emit: stats
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    seurat_qc.R \\
        $matrix_dir \\
        $prefix \\
        ${params.min_genes} \\
        ${params.min_cells} \\
        ${params.max_mito_pct} \\
        ${params.min_counts} \\
        ${params.max_counts ?: 'null'} \\
        ${params.max_genes ?: 'null'}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(as.character(getRversion()))")
        seurat: \$(Rscript -e "cat(as.character(packageVersion('Seurat')))")
    END_VERSIONS
    """
}
