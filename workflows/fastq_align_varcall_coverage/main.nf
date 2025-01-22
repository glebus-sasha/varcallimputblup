include { QC_TRIM                       } from '../qc_trim'
include { ALIGN_VARCALL                 } from '../align_varcall'
include { COVERAGE_SUMMARY              } from '../coverage_summary'
include { COVERAGE_SUMMARY_MULTIQC      } from '../../modules/multiqc/coverage_summary_multiqc'

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
    COVERAGE_SUMMARY_MULTIQC(
        QC_TRIM.out.fastp.collect(),
        QC_TRIM.out.fastqc_before.collect(),
        QC_TRIM.out.fastqc_after.collect(),
        ALIGN_VARCALL.out.flagstat.collect(),
        ALIGN_VARCALL.out.bcfstats1.map{it -> it[1]}.collect(),
        ALIGN_VARCALL.out.mosdepth.map{it -> it[1]}.collect()
    )
    emit:
    bcf = ALIGN_VARCALL.out.bcf
}