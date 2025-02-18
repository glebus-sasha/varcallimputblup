// Include processes
include { FASTQ_QC_TRIM_FASTQ_FASTP             } from '../fastq_qc_trim_fastqc_fastp'
include { BAM_VARCALL_BCFTOOLS                  } from '../bam_varcall_bcftools'
include { FASTQ_ALIGN_BWA                       } from '../fastq_align_bwa'
include { BAM_IMPUTE_GLIMPSE2                   } from '../bam_impute_glimpse2'
include { BAM_DOWNSAMPLE_SAMTOOLS               } from '../bam_downsample_samtools'
include { BCFTOOLS_INDEX                        } from '../../modules/bcftools/index'
include { BCFTOOLS_STATS                        } from '../../modules/bcftools/stats'
include { GLIMPSE2_CONCORDANCE                  } from '../../modules/glimpse2/concordance'
include { BCF_RENAME_BCFTOOLS                   } from '../../modules/local/bcf_rename_bcftools'
include { VARCALL_REF_BCFTOOLS                  } from '../../modules/local/varcall_ref_bcftools'
include { IMPUTATION_ACCURACY_PLOT              } from '../../modules/local/imputation_accuracy_plot'
include { MULTIQC                               } from '../../modules/multiqc'

workflow IMPUTATION{
    take:
    reference
    input_fastqs
    bwaidx
    faidx
    ref_panel

    main:
    FASTQ_QC_TRIM_FASTQ_FASTP(input_fastqs)
    FASTQ_ALIGN_BWA(reference, FASTQ_QC_TRIM_FASTQ_FASTP.out.trimmed_reads, bwaidx)
    BAM_DOWNSAMPLE_SAMTOOLS(FASTQ_ALIGN_BWA.out.align.map {[it[0], it[1]]})
    BAM_IMPUTE_GLIMPSE2(ref_panel, BAM_DOWNSAMPLE_SAMTOOLS.out.align)
    BCFTOOLS_INDEX(BAM_IMPUTE_GLIMPSE2.out.bcf)
    VARCALL_REF_BCFTOOLS(reference, faidx, FASTQ_ALIGN_BWA.out.align.combine(ref_panel))

    ch_validate = VARCALL_REF_BCFTOOLS.out.bcf
    ch_imputed  = BAM_IMPUTE_GLIMPSE2.out.bcf.join(BCFTOOLS_INDEX.out.csi)
    
    BCF_RENAME_BCFTOOLS(ch_validate.join(ch_imputed))
   
    ch_chrs                 = ref_panel.map{it[0]}.collect()
    ch_sids                 = ch_validate.map{it[0]}.collect()
    ch_ref_panel_with_index = ref_panel.map{[it[1], it[2]]}.collect()
    ch_validate_flat        = ch_validate
        .map{[it[0], it[2], it[3]]}
        .groupTuple()
        .map { it.flatten() }
        .map { [it[0], it[1..-1]] }
    ch_concordance = BCF_RENAME_BCFTOOLS.out.imputed_fixed_bcf.combine(ch_validate_flat, by: 0)

    GLIMPSE2_CONCORDANCE(ch_chrs, ch_ref_panel_with_index, ch_concordance)
    IMPUTATION_ACCURACY_PLOT(GLIMPSE2_CONCORDANCE.out.rsquare_grp)

    MULTIQC(
        FASTQ_QC_TRIM_FASTQ_FASTP.out.fastp_json                       |
        mix(FASTQ_QC_TRIM_FASTQ_FASTP.out.fastqc)                      |
        mix(FASTQ_QC_TRIM_FASTQ_FASTP.out.fastqc_trimmed)              |
        mix(FASTQ_ALIGN_BWA.out.flagstat.map{it[1]})                   |
       // mix(BCFTOOLS_STATS.out.bcfstats.map{it[1]})                   |
        mix(BAM_IMPUTE_GLIMPSE2.out.bcfstats_imputed.map{it[1]})       |
        collect,
        'summary'
    ) 
}