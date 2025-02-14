// Include processes
include { QC_TRIM                                   } from '../qc_trim'
include { ALIGN_VARCALL                             } from '../align_varcall'
include { BAM_VARCALL_BCFTOOLS                      } from '../bam_varcall_bcftools'
include { FASTQ_ALIGN_BWA                           } from '../fastq_align_bwa'
include { IMPUTE                                    } from '../impute'
include { BAM_DOWNSAMPLE_SAMTOOLS                   } from '../bam_downsample_samtools'
include { SAMTOOLS_VIEW as SAMTOOLS_DOWNSAMPLE      } from '../../modules/samtools/view'
include { SAMTOOLS_INDEX                            } from '../../modules/samtools/index'
include { MULTIQC                                   } from '../../modules/multiqc'

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
    
    /*ALIGN_VARCALL(
        reference,
        QC_TRIM.out.trimmed_reads,
        bwaidx,
        faidx
    )*/
    FASTQ_ALIGN_BWA(reference, QC_TRIM.out.trimmed_reads, bwaidx)
    BAM_VARCALL_BCFTOOLS(reference, faidx, FASTQ_ALIGN_BWA.out.align)

    BAM_DOWNSAMPLE_SAMTOOLS(FASTQ_ALIGN_BWA.out.align.map { it -> [it[0], it[1]] })

    IMPUTE(
        ref_panel_with_index, 
        ref_panel_index, 
        BAM_DOWNSAMPLE_SAMTOOLS.out.align
    )
    MULTIQC(
        QC_TRIM.out.fastp_json                                  |
        mix(QC_TRIM.out.fastqc)                                 |
        mix(QC_TRIM.out.fastqc_trimmed)                         |
        mix(FASTQ_ALIGN_BWA.out.flagstat.map{it -> it[1]})      |
        mix(BAM_VARCALL_BCFTOOLS.out.bcfstats.map{it -> it[1]}) |
        mix(IMPUTE.out.bcfstats_imputed.map{it -> it[1]})       |
        collect,
        'summary'
    ) 
}