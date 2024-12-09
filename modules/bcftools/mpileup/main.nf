// Define the `BCFTOOLS_MPILEUP` process that performs variant calling
process BCFTOOLS_MPILEUP {
    container = 'staphb/bcftools:latest'
    tag "$reference $bamFile ${bedfile}"
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/BCFTOOLS_MPILEUP"
//	debug true
  errorStrategy 'ignore'
	
    input:
    path reference
    tuple val(sid), path(bai), path(bamFile)
    path fai
    path bedfile

    output:
    tuple val("${sid}"), path("${sid}.vcf"),         emit: vcf
    
    script:
    def bed_option = bedfile.getBaseName() == 'dummy' ? "" : "-R  ${bedfile}"    // If the base name of bedfile is 'dummy', set bed_option to an empty string
    """    
    bcftools mpileup ${bed_option} -f $reference $bamFile -Ou | bcftools call -m -Ov -o ${sid}.vcf --threads ${task.cpus}
    """
}
