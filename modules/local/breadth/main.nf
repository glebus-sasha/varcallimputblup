// Define the `BAM_BREADTH` process that calculates the breadth of coverage
process BAM_BREADTH {
    container ''
    conda 'bioconda::bedtools'
    tag {
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
    cpus 10
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/BAM_BREADTH"
//	  debug true
//    errorStrategy 'ignore'

    input:
    tuple val(sid), path(bamFile), path(bamIndex)

    output:
    tuple val(sid), path("${bamFile.baseName}_breadth.txt"), emit: breadth

    script:
    """
    bedtools bamtobed -i ${bamFile} | \
    bedtools merge -i - | \
    awk '{sum += \$3 - \$2} END {print sum}' > ${bamFile.baseName}_breadth.txt
    """
}