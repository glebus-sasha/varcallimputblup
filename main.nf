#!/usr/bin/env nextflow

// Include workflows and modules
include { QC_TRIM       } from './workflows/qc_trim'
include { ALIGN_VARCALL } from './workflows/align_varcall'
include { IMPUTE        } from './workflows/impute'
include { MULTIQC       } from './modules/multiqc'
include { BAM_BREADTH   } from './modules/local/breadth'
include { BAM_DEPTH     } from './modules/local/depth'
include { COV_STATS     } from './modules/local/cov_stats'
include { COV_SUMMARY   } from './modules/local/cov_summary'

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
bam = Channel.fromPath("${params.bam}/*.bam").map{file->[file.simpleName, file]}
bamindex = Channel.fromPath("${params.bam}/*.bam.bai").map{file->[file.simpleName, file]}
align = bam.join(bamindex)

workflow FASTQ_ALIGN_VARCALL_COVERAGE{
    take:
    reference
    input_fastqs
    bwaidx
    faidx

    main:
    QC_TRIM(
        input_fastqs
    )
    ALIGN_VARCALL(
        reference,
        QC_TRIM.out.trimmed_reads,
        bwaidx,
        faidx
    )

    ALIGN_VARCALL.out.align |
    BAM_BREADTH & BAM_DEPTH

    breadth = BAM_BREADTH.out.breadth
    depth = BAM_DEPTH.out.depth_stats
    bcfstats1 = ALIGN_VARCALL.out.bcfstats1

    COV_STATS(breadth.join(depth).join(bcfstats1))
    COV_SUMMARY(COV_STATS.out.cov_stats.map{it -> it[1]}.collect())

    
}

workflow imputation{
    take:
    reference
    input_fastqs
    bwaidx
    faidx
    ref_panel_with_index
    ref_panel_index


    main:
    QC_TRIM(
        input_fastqs
    )
    ALIGN_VARCALL(
        reference,
        QC_TRIM.out.trimmed_reads,
        bwaidx,
        faidx
    )
    IMPUTE(
        ref_panel_with_index, 
        ref_panel_index, 
        ALIGN_VARCALL.out.align
    )
    MULTIQC(
        QC_TRIM.out.fastp.collect(),
        QC_TRIM.out.fastqc_before.collect(),
        QC_TRIM.out.fastqc_after.collect(),
        ALIGN_VARCALL.out.flagstat.collect(),
        ALIGN_VARCALL.out.bcfstats1.collect(),
        IMPUTE.out.bcfstats2.collect()
    )
}

workflow {
    FASTQ_ALIGN_VARCALL_COVERAGE(
        reference,
        input_fastqs,
        bwaidx,
        faidx
    )
}

