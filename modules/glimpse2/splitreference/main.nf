process GLIMPSE2_SPLITREFERENCE {
    tag "${ref_panel.baseName}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container 'imary116/glimpse2:with-bcftools-and-updated-info-score'
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/GLIMPSE2_SPLITREFERENCE"

    input:
    path(ref_panel)
    path(ref_panel_index)
    val(chunk_chr)
    val(input_region)
    val(output_region)


    output:
        tuple val(meta), path("*.bin"), emit: bin_ref
        path "versions.yml"           , emit: versions

    script:
    def prefix      = task.ext.prefix ?: "${ref_panel.baseName}_${output_region.replace(":","_")}"

    shell:
    """
    REF=!ref_panel
    while IFS="" read -r LINE || [ -n "!LINE" ];
    do
        printf -v ID "%02d" $(echo !LINE | cut -d" " -f1)
        IRG=$(echo !LINE | cut -d" " -f3)
        ORG=$(echo !LINE | cut -d" " -f4)

        ./bin/GLIMPSE2_split_reference --reference !{REF} --input-region !{IRG} --output-region !{ORG} --output reference_panel/split/1000GP.chr22.noNA12878
    done < !chunk_chr

    stub:
    """
        GLIMPSE2_split_reference \
        --reference $ref_panel \
        --input-region $input_region \
        --output-region $output_region \
        --thread $task.cpus \
        --output ${prefix}
    touch ${ref_panel.baseName}.bin
    """
}
