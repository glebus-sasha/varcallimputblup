#!/usr/bin/env nextflow

// Include workflows
include { CLUSTERING } from './workflows/clustering'
include { IMPUTATION } from './workflows/imputation'


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
ref_panel = Channel.fromPath("${params.ref_panel}").map{file->[file.simpleName, file]}
ref_panel_index = Channel.fromPath("${params.ref_panel_index}").map{file->[file.simpleName, file]}
ref_panel_with_index = ref_panel.join(ref_panel_index)

// aligments channels
/* bam = Channel.fromPath("${params.bam}/*.bam").map{file->[file.simpleName, file]}
bamindex = Channel.fromPath("${params.bam}/*.bam.bai").map{file->[file.simpleName, file]}
align = bam.join(bamindex)
*/

workflow imputation{
    take:
    reference
    input_fastqs
    bwaidx
    faidx
    ref_panel_with_index
    ref_panel_index

    main:
    IMPUTATION(
        reference,
        input_fastqs,
        bwaidx,
        faidx,
        ref_panel_with_index,
        ref_panel_index
    )
}

workflow clustering{
    take:
    reference
    input_fastqs
    bwaidx
    faidx

    main:
    CLUSTERING(
        reference,
        input_fastqs,
        bwaidx,
        faidx
    )
}

workflow{
    clustering(
        reference,
        input_fastqs,
        bwaidx,
        faidx
    )
}