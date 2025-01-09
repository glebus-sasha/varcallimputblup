// Define the `BAM_DEPTH` process that calculates the depth of coverage
process BAM_DEPTH {
    container ''
    conda 'bioconda::samtools'
    tag {
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
    cpus 10
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/BAM_DEPTH"
//	  debug true
//    errorStrategy 'ignore'

    input:
    tuple val(sid), path(bamFile), path(bamIndex)

    output:
    tuple val(sid), path("${bamFile.baseName}_depth.txt"), emit: depth

    script:
    """
    samtools depth -a ${bamFile} > ${bamFile.baseName}_depth.txt
    """
}
