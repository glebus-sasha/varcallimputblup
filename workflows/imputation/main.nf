// Include processes
include { FASTQ_QC_TRIM_FASTQ_FASTP             } from '../fastq_qc_trim_fastqc_fastp'
include { FASTQ_ALIGN_BWA                       } from '../fastq_align_bwa'
include { BAM_IMPUTE_GLIMPSE2                   } from '../bam_impute_glimpse2'
include { BAM_DOWNSAMPLE_SAMTOOLS               } from '../bam_downsample_samtools'
include { BCF_ACCURACY_GLIMPSE2                 } from '../bcf_accuracy_glimpse2'
include { MOSDEPTH                              } from '../../modules/mosdepth'
include { MULTIQC                               } from '../../modules/multiqc'

workflow IMPUTATION{
    take:
    reference
    input_fastqs
    bwaidx
    faidx
    ref_panel
    downsample_rate

    main:
    FASTQ_QC_TRIM_FASTQ_FASTP(input_fastqs)
    FASTQ_ALIGN_BWA(reference, FASTQ_QC_TRIM_FASTQ_FASTP.out.trimmed_reads, bwaidx)
    BAM_DOWNSAMPLE_SAMTOOLS(FASTQ_ALIGN_BWA.out.align.map {[it[0], it[1]]}, downsample_rate)
    BAM_IMPUTE_GLIMPSE2(ref_panel, BAM_DOWNSAMPLE_SAMTOOLS.out.align)
    MOSDEPTH(BAM_DOWNSAMPLE_SAMTOOLS.out.align, 'downsampled')
    BCF_ACCURACY_GLIMPSE2(
        reference, 
        faidx, 
        BAM_IMPUTE_GLIMPSE2.out.bcf, 
        FASTQ_ALIGN_BWA.out.align, 
        ref_panel
        )

    MULTIQC(
        FASTQ_QC_TRIM_FASTQ_FASTP.out.fastp_json                       |
        mix(FASTQ_QC_TRIM_FASTQ_FASTP.out.fastqc)                      |
        mix(FASTQ_QC_TRIM_FASTQ_FASTP.out.fastqc_trimmed)              |
        mix(FASTQ_ALIGN_BWA.out.flagstat.map{it[1]})                   |
        mix(FASTQ_ALIGN_BWA.out.mosdepth.map{it[1]})                   |
        mix(FASTQ_ALIGN_BWA.out.mosdepth_summary.map{it[1]})           |
        mix(BAM_IMPUTE_GLIMPSE2.out.bcfstats_imputed.map{it[1]})       |
        mix(BCF_ACCURACY_GLIMPSE2.out.glimpse2_errors_grp.map{it[1]})  |
        mix(BCF_ACCURACY_GLIMPSE2.out.glimpse2_errors_spl.map{it[1]})  |
        mix(MOSDEPTH.out.global_dist.map{it[1]})                       |
        mix(MOSDEPTH.out.summary.map{it[1]})                           |
        collect,
        'summary'
    ) 
}