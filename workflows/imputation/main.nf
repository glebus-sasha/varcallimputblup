// Include processes
include { QC_TRIM       } from '../qc_trim'
include { ALIGN_VARCALL } from '../align_varcall'
include { IMPUTE        } from '../impute'
include { MULTIQC       } from '../../modules/multiqc'

workflow IMPUTATION{
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
        QC_TRIM.out.fastp                                   |
        mix(QC_TRIM.out.fastqc_before)                      |
        mix(QC_TRIM.out.fastqc_after)                       |
        mix(ALIGN_VARCALL.out.flagstat)                     |
        mix(ALIGN_VARCALL.out.bcfstats1.map{it -> it[1]})   |
        mix(IMPUTE.out.bcfstats2.map{it -> it[1]})          |
        collect,
        'summary'
    ) 
}