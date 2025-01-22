// Define the `COMBINE_GDS` process that performs cobining gds files
process COMBINE_GDS {
    conda "${moduleDir}/environment.yml"
    tag 'all_samples'
//    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/COMBINE_GDS"
//	debug true
    errorStrategy 'ignore'
	
    input:
    path gds_files

    output:
    path 'combined.gds', emit: combined_gds

    script:
    """    
    #!/usr/bin/env Rscript

    library(gdsfmt)
    library(SNPRelate)

    # List all .gds files in the current directory
    gds_files <- list.files(full.names = TRUE)

    # Combine the GDS files
    snpgdsCombineGeno(gds_files, 'combined.gds')
    """
}