// Define the `COV_SUMMARY` process that generates a summary table for all samples
process COV_SUMMARY {
    container ''
    conda "${moduleDir}/environment.yml"
    tag 'summary'
    cpus 10
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/COV_SUMMARY"
//	  debug true
//    errorStrategy 'ignore'

    input:
    tuple val(sid), path(breadthFile), path(depthStatsFile), path(bcfstatsFile)

    output:
    path("combined_summary.txt")

    script:
    """
    #!/usr/bin/env Rscript

    library(dplyr)
    library(purrr)
    library(readr)

    # Function to read and process a single sample's stats files
    process_sample <- function(sid, depthStatsFile, breadthFile, bcfstatsFile) {
        depth_stats <- read_tsv(depthStatsFile)
        breadth <- read_tsv(breadthFile)
        bcfstats <- read_tsv(bcfstatsFile)

        # Extract the number of SNPs
        snp_count <- as.numeric(bcfstats[bcfstats\$Type == "SNP", "Count"])

        # Combine the stats into a single data frame
        sample_stats <- tibble(
            Sample = sid,
            MeanDepth = depth_stats\$Mean,
            MedianDepth = depth_stats\$Median,
            MinDepth = depth_stats\$Min,
            MaxDepth = depth_stats\$Max,
            SDDepth = depth_stats\$SD,
            Breadth = breadth\$V1,
            SNPCount = snp_count
        )

        return(sample_stats)
    }

    # Combine all samples' stats into a single data frame
    combined_results <- pmap_dfr(
        list(${sid}, ${depthStatsFile}, ${breadthFile}, ${bcfstatsFile}),
        process_sample
    )

    # Write the combined results to a file
    write_tsv(combined_results, "combined_summary.txt")
    """
}