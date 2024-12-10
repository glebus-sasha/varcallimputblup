process GLIMPSE2_CHUNK {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/glimpse-bio:2.0.1--h46b9e50_1':
        'biocontainers/glimpse-bio:2.0.1--h46b9e50_1' }"
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/GLIMPSE2_CHUNK"

    input:
    path(ref_panel)
    path(ref_panel_index)
    path(region)

    output:
    path("${ref_panel.baseName}.txt"), emit: chunk_chr

    script:
    """
    GLIMPSE2_chunk \
        --input $input \
        --region $region \
        --threads $task.cpus \
        --output ${ref_panel.baseName}.txt
    """

    stub:
    """
    touch ${ref_panel.baseName}.txt
    """
}
