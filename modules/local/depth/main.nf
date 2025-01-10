// Define the `BAM_DEPTH` process that calculates the depth of coverage and generates a plot and summary statistics
process BAM_DEPTH {
    container 'glebusasha/r_env_image:latest'
    conda "${moduleDir}/environment.yml"
    tag {
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/BAM_DEPTH", pattern: '*.png'
//	  debug true
    errorStrategy 'ignore'

    input:
    tuple val(sid), path(bamFile), path(bamIndex)

    output:
    tuple val(sid), path("${sid}_depth_plot.png") , emit: depth_plot
    tuple val(sid), path("${sid}_depth_stats.csv"), emit: depth_stats

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

    p <- ggplot(bam, aes(x=pos, y=depth)) +
    geom_line() +
    labs(title="Depth of Coverage", x="Genomic Position", y="Depth") +
    theme_bw() +
    scale_x_continuous(breaks = seq(0, 1e9, by = 1e7), labels = scales::comma) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave("${sid}_depth_plot.png", p, width=10, height=6)

    summary_stats <- bam %>%
      summarise(
        Mean_Depth = mean(depth, na.rm = TRUE),
        Median_Depth = median(depth, na.rm = TRUE),
        Min_Depth = min(depth, na.rm = TRUE),
        Max_Depth = max(depth, na.rm = TRUE),
        SD_Depth = sd(depth, na.rm = TRUE)
      )
    write.csv(summary_stats, file="${sid}_depth_stats.csv", row.names=FALSE, quote=FALSE)
    """
}
