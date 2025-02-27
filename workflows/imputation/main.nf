// Include processes
include { FASTQ_QC_TRIM_FASTQ_FASTP             } from '../../subworkflows/fastq_qc_trim_fastqc_fastp'
include { FASTQ_ALIGN_BWA                       } from '../../subworkflows/fastq_align_bwa'
include { BAM_IMPUTE_GLIMPSE2                   } from '../../subworkflows/bam_impute_glimpse2'
include { BAM_DOWNSAMPLE_SAMTOOLS               } from '../../subworkflows/bam_downsample_samtools'
include { BCF_ACCURACY_GLIMPSE2                 } from '../../subworkflows/bcf_accuracy_glimpse2'
include { CARPET                                } from '../../subworkflows/carpet'
include { MOSDEPTH                              } from '../../modules/mosdepth'
include { PICARD_MARK_DUPLICATES                } from '../../modules/picard/mark_duplicates'
include { MULTIQC                               } from '../../modules/multiqc'

workflow IMPUTATION{
    take:
    reference
    input_fastqs
    bwaidx
    faidx
    ref_panel
    bed

    main:
    FASTQ_QC_TRIM_FASTQ_FASTP(input_fastqs)
    FASTQ_ALIGN_BWA(reference, FASTQ_QC_TRIM_FASTQ_FASTP.out.trimmed_reads, bwaidx)
    PICARD_MARK_DUPLICATES(FASTQ_ALIGN_BWA.out.align)
    BAM_IMPUTE_GLIMPSE2(ref_panel, FASTQ_ALIGN_BWA.out.align)
    CARPET(BAM_IMPUTE_GLIMPSE2.out.bcf.join(BAM_IMPUTE_GLIMPSE2.out.csi), bed)

    MULTIQC(
        FASTQ_QC_TRIM_FASTQ_FASTP.out.fastp_json                       |
        mix(FASTQ_QC_TRIM_FASTQ_FASTP.out.fastqc)                      |
        mix(FASTQ_QC_TRIM_FASTQ_FASTP.out.fastqc_trimmed)              |
        mix(FASTQ_ALIGN_BWA.out.flagstat.map{it[1]})                   |
        mix(FASTQ_ALIGN_BWA.out.mosdepth.map{it[1]})                   |
        mix(FASTQ_ALIGN_BWA.out.mosdepth_summary.map{it[1]})           |
        mix(BAM_IMPUTE_GLIMPSE2.out.bcfstats_imputed.map{it[1]})       |
        mix(PICARD_MARK_DUPLICATES.out.metrics.map{it[1]})             |
        collect,
        'summary'
    ) 
}