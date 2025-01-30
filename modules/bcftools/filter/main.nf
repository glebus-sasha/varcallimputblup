// Define the `BCFTOOLS_FILTER` process that performs vcf statistics
process BCFTOOLS_FILTER {
    container 'staphb/bcftools:latest'
    conda 'bcftools'
    tag { 
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
//    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/BCFTOOLS_FILTER"

//    cache "lenient" 
//    debug true
    errorStrategy 'ignore'
	
    input:
    tuple val(sid), path(bcf)

    output:
    tuple val(sid), path("${sid}_filtered.bcf"),      emit: bcf
    
    script:
    """
    bcftools filter -i 'INFO/DP > 2' ${bcf} -o "${sid}_filtered.bcf"
    """
}