// Define the `GDS_CLUSTER` process that performs clustering
process GDS_CLUSTER {
    conda "${moduleDir}/environment.yml"
    tag 'all_samples'
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/GDS_CLUSTER"
//	debug true
    errorStrategy 'ignore'
	
    input:
    path gds_file

    output:
    path 'cluster_tab.tsv'  , emit: cluster_tab
    path 'dendrogram.nwk'   , emit: dendrogram_nwk
    path 'dendrogram.png'   , emit: dendrogram_png
    path 'pca.png'          , emit: pca_png
    path 'pca_plot.png'     , emit: pca_plot_png

    script:
    """    
    #!/usr/bin/env Rscript

    source('${projectDir}/assets/cluster.R')
    gds_cluster()
    """
}