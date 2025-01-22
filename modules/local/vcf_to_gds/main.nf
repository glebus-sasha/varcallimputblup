// Define the `VCF_TO_GDS` process that converts BCF files to GDS format
process VCF_TO_GDS {
    conda "${moduleDir}/environment.yml"
    tag { 
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
//    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/VCF_TO_GDS"
//    debug true
    errorStrategy 'ignore'

    input:
    tuple val(sid), path(vcf)

    output:
    tuple val(sid), path("${sid}.gds"), emit: gds

    script:
    """
    #!/usr/bin/env Rscript

    library(SNPRelate)
    library(SeqArray)
    snpgdsVCF2GDS("${vcf}", "${sid}.gds", method = "biallelic.only")
    """
}
