include { QC_TRIM          } from '../qc_trim'
include { ALIGN_VARCALL    } from '../align_varcall'
include { COVERAGE_SUMMARY } from '../coverage_summary'
include { BCF_CLUSTERING   } from '../bcf_clustering'
include { MULTIQC          } from '../../modules/multiqc'

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
        ALIGN_VARCALL.out.mosdepth_dedup_summary, 
        reference
    )
    BCF_CLUSTERING(ALIGN_VARCALL.out.bcf)
    MULTIQC(
        QC_TRIM.out.fastp_json                                  |
        mix(QC_TRIM.out.fastqc)                                 |
        mix(QC_TRIM.out.fastqc_trimmed)                         |
        mix(ALIGN_VARCALL.out.flagstat.map{it -> it[1]})        |
        mix(ALIGN_VARCALL.out.flagstat_dedup.map{it -> it[1]})  |
        mix(ALIGN_VARCALL.out.bcfstats.map{it -> it[1]})        |
        mix(ALIGN_VARCALL.out.mosdepth.map{it -> it[1]})        |
        mix(ALIGN_VARCALL.out.mosdepth_dedup.map{it -> it[1]})  |
        collect,
        'summary'
    )
}