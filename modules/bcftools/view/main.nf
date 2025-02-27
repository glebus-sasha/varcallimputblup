// Define the `BCFTOOLS_VIEW` process that performs vcf statistics
process BCFTOOLS_VIEW {
    container 'staphb/bcftools:latest'
    conda 'bcftools'
    tag { 
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/BCFTOOLS_VIEW"

//    cache "lenient" 
//    debug true
//    errorStrategy 'ignore'
	
    input:
    tuple val(sid), path(bcf), path(csi)
    path bed
    val tag

    output:
    tuple val(sid), path("${sid}${tag}.vcf"),      emit: vcf
    
    script:
    """
    bcftools view -R ${bed} ${bcf} >  "${sid}${tag}.vcf"
    """
}