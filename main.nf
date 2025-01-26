#!/usr/bin/env nextflow

// Include workflows and modules
include { QC_TRIM                       } from './workflows/qc_trim'
include { ALIGN_VARCALL                 } from './workflows/align_varcall'
include { IMPUTE                        } from './workflows/impute'
include { COVERAGE_SUMMARY              } from './workflows/coverage_summary'
include { CUTADAPT_QC                   } from './workflows/cutadapt_qc'
include { BCF_CLUSTERING                } from './workflows/bcf_clustering'
include { FASTQ_ALIGN_VARCALL_COVERAGE  } from './workflows/fastq_align_varcall_coverage'
include { MULTIQC                       } from './modules/multiqc'

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
    QC_TRIM.out.fastp                                       |
        mix(QC_TRIM.out.fastqc_before)                      |
        mix(QC_TRIM.out.fastqc_after)                       |
        mix(ALIGN_VARCALL.out.flagstat)                     |
        mix(ALIGN_VARCALL.out.bcfstats1.map{it -> it[1]})   |
        mix(IMPUTE.out.bcfstats2.map{it -> it[1]})          |
        collect                                             |
        MULTIQC
}

workflow {
    FASTQ_ALIGN_VARCALL_COVERAGE(
        reference,
        input_fastqs,
        bwaidx,
        faidx
    )
    BCF_CLUSTERING(FASTQ_ALIGN_VARCALL_COVERAGE.out.bcf)
}