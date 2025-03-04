// Define the `BCF_TO_VCF` process that converts BCF files to GDS format
process BCF_TO_VCF {
    conda "bcftools"
    tag { 
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
//    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/BCF_TO_VCF"
//    debug true
    errorStrategy 'ignore'

    input:
    tuple val(sid), path(bcf), path(csi)

    output:
    tuple val(sid), path("${sid}.vcf"), emit: vcf

    script:
    """
    bcftools view -Ov ${bcf} -o ${sid}.vcf
    """
}
