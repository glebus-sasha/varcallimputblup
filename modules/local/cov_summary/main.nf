// Define the `COV_SUMMARY` process that generates a summary coverage table for each sample
process COV_SUMMARY {
    container ''
    conda 'r-base'
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
    breadth <- as.numeric(readLines("${breadthFile}"))
    depth_data <- read.table("${depthFile}", header=FALSE)
    average_depth <- mean(depth_data[,3])
    summary_data <- data.frame(Sample=c("${sid}"), Breadth=breadth, AverageDepth=average_depth)
    write.table(summary_data, file="${sid}_summary.txt", sep="\\t", row.names=FALSE, col.names=TRUE)
    """
}