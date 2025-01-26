// Define the `SAMTOOLS_FLAGSTAT` process that aligns stats
process SAMTOOLS_FLAGSTAT {
    container 'glebusasha/bwa_samtools'
    conda 'bioconda::bwa bioconda::samtools'
    tag {
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
//    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/SAMTOOLS_FLAGSTAT"
//	  debug true
//    errorStrategy 'ignore'

    input:
    tuple val(sid), path(bamFile)
    val tag

    
    output:
    path "${sid}${tag}.flagstat", emit: flagstat
    
    script:
    """
    samtools flagstat $bamFile > ${sid}${tag}.flagstat
    """
}
