process VCFLIB_VCFFIXUP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcflib:1.0.3--hecb563c_1':
        'biocontainers/vcflib:1.0.3--hecb563c_1' }"

    input:
    tuple val(chr), path(ref_panel), path(ref_panel_index)

    output:
    tuple val(chr), path("${chr}_fixed.vcf.gz"), emit: vcf_gz

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    vcffixup \\
        ref_panel | bgzip -c > ${chr}_fixed.vcf.gz
    """
}