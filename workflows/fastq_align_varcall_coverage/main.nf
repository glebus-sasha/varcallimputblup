include { QC_TRIM                       } from '../qc_trim'
include { ALIGN_VARCALL                 } from '../align_varcall'
include { COVERAGE_SUMMARY              } from '../coverage_summary'
include { MULTIQC                       } from '../../modules/multiqc'

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
    COVERAGE_SUMMARY(ALIGN_VARCALL.out.align, ALIGN_VARCALL.out.bcfstats1, ALIGN_VARCALL.out.mosdepth_summary, reference)
    MULTIQC(
        QC_TRIM.out.fastp                                   |
        mix(QC_TRIM.out.fastqc_before)                      |
        mix(QC_TRIM.out.fastqc_after)                       |
        mix(ALIGN_VARCALL.out.flagstat)                     |
        mix(ALIGN_VARCALL.out.bcfstats1.map{it -> it[1]})   |
        mix(ALIGN_VARCALL.out.mosdepth.map{it -> it[1]})    |
        collect
    )
    emit:
    bcf = ALIGN_VARCALL.out.bcf
}