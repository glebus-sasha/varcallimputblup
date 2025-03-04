// Define the `VCF_CARPET` process that merges all vcf files
process VCF_CARPET {
    container 'glebusasha/r_env_image:latest'
    conda "${moduleDir}/environment.yml"
    tag 'all samples'
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/VCF_CARPET/${tag}"
    errorStrategy 'ignore'

    input:
    path vcf_files
    path bed
    path interpretation
    val tag

    output:
    path "carpet.html", emit: carpet
    path "table.csv"  , emit: table_csv
    path "table.html" , emit: table_html

    script:
    """
    Rscript ${projectDir}/assets/carpet.R \
        --bed_file ${bed} \
        --interpretation_file ${interpretation}
    """
}