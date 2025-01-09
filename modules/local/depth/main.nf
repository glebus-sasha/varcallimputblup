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
    errorStrategy 'ignore'

    input:
    tuple val(sid), path(bamFile), path(bamIndex)

    output:
    tuple val(sid), path("${sid}_depth_plot.png"), path("${sid}_depth_stats.txt"), emit: depth

    script:
    """
    #!/usr/bin/env Rscript
    library(Rsamtools)
    library(GenomicRanges)
    library(ggplot2)
    library(dplyr)

    bamFile <- "${bamFile}"
    bam <- BamFile(bamFile, asMates=TRUE) %>%
      scanBam() %>%
      as.data.frame() %>%
      group_by(rname, pos) %>%
      summarise(depth = n(), .groups = 'drop')

    p <- bam %>% ggplot() +
      aes(x=pos, y=depth) +
      geom_line() +
      labs(title="Depth of Coverage", x="Genomic Position", y="Depth") +
      theme_minimal()
    ggsave("${sid}_depth_plot.png", p, width=10, height=6)

    summary_stats <- bam %>%
      summarise(
        Mean = mean(depth, na.rm = TRUE),
        Median = median(depth, na.rm = TRUE),
        Min = min(depth, na.rm = TRUE),
        Max = max(depth, na.rm = TRUE),
        SD = sd(depth, na.rm = TRUE)
      )
    write.table(summary_stats, file="${sid}_depth_stats.txt", sep="\\t", row.names=FALSE, col.names=TRUE)
    """
}
