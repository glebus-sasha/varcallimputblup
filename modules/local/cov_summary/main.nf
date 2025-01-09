// Define the `COV_SUMMARY` process that generates a summary table for all samples
process COV_SUMMARY {
    container ''
    conda 'environment.yml'
    tag 'summary'
    cpus 10
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/COV_SUMMARY"
//	  debug true
//    errorStrategy 'ignore'

    input:
    tuple val(sid), path(depthStatsFile), path(breadthFile), path(bcfstatsFile)

    output:
    path("combined_summary.txt")

    script:
    """
    library(dplyr)
    library(readr)

    # Initialize an empty data frame to store combined results
    combined_results <- data.frame()

    # Loop through each sample and read the stats files
    for (i in 1:length(sid)) {
        depth_stats <- read_tsv(depthStatsFile[i])
        breadth <- read_tsv(breadthFile[i])
        bcfstats <- read_tsv(bcfstatsFile[i])

        # Extract the number of SNPs
        snp_count <- as.numeric(bcfstats[bcfstats\$Type == "SNP", "Count"])

        # Combine the stats into a single data frame
        sample_stats <- data.frame(
            Sample = sid[i],
            MeanDepth = depth_stats\$Mean,
            MedianDepth = depth_stats\$Median,
            MinDepth = depth_stats\$Min,
            MaxDepth = depth_stats\$Max,
            SDDepth = depth_stats\$SD,
            Breadth = breadth\$V1,
            SNPCount = snp_count
        )

        # Append to the combined results
        combined_results <- bind_rows(combined_results, sample_stats)
    }

    # Write the combined results to a file
    write_tsv(combined_results, "combined_summary.txt")
    """
}
