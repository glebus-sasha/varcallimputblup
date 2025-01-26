// Define the `SAMTOOLS_RMDUP` process that removes duplicates
process SAMTOOLS_RMDUP {
    container 'glebusasha/bwa_samtools'
    conda 'bioconda::bwa bioconda::samtools'
    tag {
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
   publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/SAMTOOLS_RMDUP"
//	  debug true
//    errorStrategy 'ignore'
    input:
    tuple val(sid), path(bamFile)

    output:
    tuple val(sid), path("${sid}_dedup.bam"), emit: bam

    script:
    """
    samtools rmdup ${bamFile} "${sid}_dedup.bam"
    """
}