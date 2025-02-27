#!/usr/bin/env nextflow

// Include workflows
include { CLUSTERING                     } from './workflows/clustering'
include { IMPUTATION_LOW_PASS_SIMULATION } from './workflows/imputation_low_pass_simulation'
include { IMPUTATION                     } from './workflows/imputation'


// Logging pipeline information
log.info """\
\033[0;36m  ==========================================  \033[0m
\033[0;34m       v a r c a l l i m p u t b l u p        \033[0m
\033[0;36m  ==========================================  \033[0m
    """

// reference and reads channels
reference = Channel.fromPath("${params.reference}").collect()
input_fastqs = Channel.fromFilePairs(["${params.reads}/*[rR]{1,2}*.*{fastq,fq}*", 
                                      "${params.reads}/*_{1,2}.{fastq,fq}*"], flat: true)
// index channels
bwaidx = Channel.fromPath("${params.bwaidx}/*").collect()
faidx = Channel.fromPath("${params.faidx}/*.fai").collect()

// reference panel channels
ref_panel_vcf       = Channel.fromPath("${params.ref_panel_vcf}").map{file->[file.simpleName, file]}
ref_panel_index     = Channel.fromPath("${params.ref_panel_index}").map{file->[file.simpleName, file]}
ref_panel_tsv       = Channel.fromPath("${params.ref_panel_tsv}").map{file->[file.simpleName, file]}
ref_panel_tsv_index = Channel.fromPath("${params.ref_panel_tsv_index}").map{file->[file.simpleName, file]}
ref_panel = ref_panel_vcf.join(ref_panel_index).join(ref_panel_tsv).join(ref_panel_tsv_index)

downsample_rate = params.downsample_rate
// aligments channels
/* bam = Channel.fromPath("${params.bam}/*.bam").map{file->[file.simpleName, file]}
bamindex = Channel.fromPath("${params.bam}/*.bam.bai").map{file->[file.simpleName, file]}
align = bam.join(bamindex)
*/



workflow{
    /*
    CLUSTERING(
        reference,
        input_fastqs,
        bwaidx,
        faidx
    )
    */
    /*
    IMPUTATION_LOW_PASS_SIMULATION(
        reference,
        input_fastqs,
        bwaidx,
        faidx,
        ref_panel,
        downsample_rate
    )*/


    // regions channels
    bed = Channel.fromPath("${params.bed}").collect()
    IMPUTATION(
        reference,
        input_fastqs,
        bwaidx,
        faidx,
        ref_panel,
        bed
    )
}