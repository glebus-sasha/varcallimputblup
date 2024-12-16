// Define the `BCFTOOLS_INDEX` process that performs vcf statistics
process BCFTOOLS_INDEX {
    container = 'staphb/bcftools:latest'
    tag "$sid"
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/BCFTOOLS_INDEX"

//    cache "lenient" 
//    debug true
    errorStrategy 'ignore'
	
    input:
    tuple val(sid), path(bcf)

    output:
    path "${sid}.bcf.csi",      emit: index
    
    script:
    """
    bcftools index $bcf
    """
}