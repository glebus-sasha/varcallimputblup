// Define the `BAM_DEPTH` process that calculates the depth of coverage and generates a plot and summary statistics
process BAM_DEPTH {
    container ''
    conda "${moduleDir}/environment.yml"
    tag {
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
    cpus 10
//    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/BAM_DEPTH"
//	  debug true
//    errorStrategy 'ignore'

    input:
    tuple val(sid), path(bamFile), path(bamIndex)

    output:
    tuple val(sid), path("${sid}_depth_plot.png"), path("${sid}_depth_stats.txt"), emit: depth

    script:
    """
    library(Rsamtools)
    library(GenomicRanges)
    library(ggplot2)

    bamFile <- "${bamFile}"
    bam <- BamFile(bamFile, asMates=TRUE)
    coverage <- as.data.frame(scanBam(bam, which=c("qwidth", "qstart", "qend")))

    # Calculate depth of coverage
    depth <- coverage\$qend - coverage\$qstart + 1
    coverage\$depth <- depth

    # Generate plot
    p <- ggplot(coverage, aes(x=qstart, y=depth)) +
        geom_line() +
        labs(title="Depth of Coverage", x="Genomic Position", y="Depth") +
        theme_minimal()
    ggsave("${sid}_depth_plot.png", p, width=10, height=6)

    # Calculate summary statistics
    summary_stats <- data.frame(
        Sample = "${sid}",
        MeanDepth = mean(depth),
        MedianDepth = median(depth),
        MinDepth = min(depth),
        MaxDepth = max(depth),
        SDDepth = sd(depth)
    )
    write.table(summary_stats, file="${sid}_depth_stats.txt", sep="\\t", row.names=FALSE, col.names=TRUE)
    """
}
