// Define the `SAMTOOLS_VIEW` process that performs samtools view
process SAMTOOLS_VIEW {
    container 'glebusasha/bwa_samtools'
    conda 'bioconda::bwa bioconda::samtools'
    tag {
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
//    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/SAMTOOLS_VIEW"
//	  debug true
//    errorStrategy 'ignore'

    input:
    tuple val(sid), path(bamFile)
    val freq

    
    output:
    tuple val(sid), path("${sid}_sim_${freq}.bam"), emit: bam
    
    script:
    """
    samtools view -s $freq -b $bamFile > ${sid}_sim_${freq}.bam
    """
}
