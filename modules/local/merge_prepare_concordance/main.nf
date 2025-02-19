// Define the `MERGE_PREPARE_CONCORDANCE` process that plots concordance accuracy
process MERGE_PREPARE_CONCORDANCE {
    conda ""
    tag 'all_samples'
//  publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/MERGE_PREPARE_CONCORDANCE"
//	debug true
    errorStrategy 'ignore'
	
    input:
    val sids
    path concordance

    output:
    path "sids.txt"                 , emit: sids
    path "merged_concordance.txt"   , emit: merged_concordance
    
    script:
    def sid_list = sids.collect { sid -> 
        """
        echo ${sid}.bam >> sids.txt
        """
        }.join("\n")

    def concordance_list = concordance.collect { concordance -> 
        """
        cat ${concordance} >> merged_concordance.txt
        """
        }.join("\n")

    
    """
    ${sid_list}
    ${concordance_list}
    """
}