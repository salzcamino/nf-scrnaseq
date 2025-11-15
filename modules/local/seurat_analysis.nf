process SEURAT_ANALYSIS {
    tag "$meta.id"
    label 'process_high'

    conda "conda-forge::r-base=4.3.1 bioconda::r-seurat=4.3.0 conda-forge::r-ggplot2=3.4.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-seurat:4.3.0' :
        'satijalab/seurat:4.3.0' }"

    input:
    tuple val(meta), path(filtered_rds)

    output:
    tuple val(meta), path("*_processed.rds")     , emit: processed
    tuple val(meta), path("*_markers.csv")       , emit: markers
    tuple val(meta), path("*_top_markers.csv")   , emit: top_markers
    tuple val(meta), path("*_cluster_stats.csv") , emit: cluster_stats
    tuple val(meta), path("*.pdf")               , emit: plots
    tuple val(meta), path("*_summary.txt")       , emit: summary
    tuple val(meta), path("*.h5ad")              , optional:true, emit: h5ad
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    seurat_analysis.R \\
        $filtered_rds \\
        $prefix \\
        ${params.n_hvgs} \\
        ${params.n_pcs} \\
        ${params.resolution} \\
        ${params.normalization_method}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(as.character(getRversion()))")
        seurat: \$(Rscript -e "cat(as.character(packageVersion('Seurat')))")
    END_VERSIONS
    """
}
