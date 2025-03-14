include { FASTQ_QC_TRIM_FASTQ_FASTP } from '../../subworkflows/fastq_qc_trim_fastqc_fastp'
include { ALIGN_VARCALL             } from '../../subworkflows/align_varcall'
include { COVERAGE_SUMMARY          } from '../../subworkflows/coverage_summary'
include { BCF_CLUSTERING            } from '../../subworkflows/bcf_clustering'
include { BCFTOOLS_FILTER           } from '../../modules/bcftools/filter'
include { MULTIQC                   } from '../../modules/multiqc'

workflow CLUSTERING {
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
    COVERAGE_SUMMARY(
        ALIGN_VARCALL.out.align, 
        ALIGN_VARCALL.out.bcfstats, 
        ALIGN_VARCALL.out.mosdepth_summary, 
        reference
    )
    BCFTOOLS_FILTER(ALIGN_VARCALL.out.bcf)
    BCF_CLUSTERING(BCFTOOLS_FILTER.out.bcf)
    MULTIQC(
        QC_TRIM.out.fastp_json                                  |
        mix(QC_TRIM.out.fastqc)                                 |
        mix(QC_TRIM.out.fastqc_trimmed)                         |
        mix(ALIGN_VARCALL.out.flagstat.map{it -> it[1]})        |
        mix(ALIGN_VARCALL.out.bcfstats.map{it -> it[1]})        |
        mix(ALIGN_VARCALL.out.mosdepth.map{it -> it[1]})        |
        collect,
        'summary'
    )
}