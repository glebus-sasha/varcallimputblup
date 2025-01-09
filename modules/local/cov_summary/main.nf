// Define the `COV_SUMMARY` process that generates a summary coverage table for each sample
process COV_SUMMARY {
    container ''
    conda ''
    tag {
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
    cpus 10
//    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/COV_SUMMARY"
//	  debug true
//    errorStrategy 'ignore'

    input:
    tuple val(sid), path(breadthFile), path(depthFile)

    output:
    tuple val(sid), path("${sid}_summary.txt"), emit: summary

    script:
    """
    breadth=\$(cat ${breadthFile})
    depth=\$(awk '{sum += \$3} END {if (NR > 0) print sum/NR; else print 0}' ${depthFile})
    echo -e "Sample\\tBreadth\\tAverageDepth\\n${sid}\\t${breadth}\\t${depth}" > ${sid}_summary.txt
    """
}