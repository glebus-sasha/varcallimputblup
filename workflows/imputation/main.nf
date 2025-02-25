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
include { PREPARE_CONCORDANCE                   } from '../../modules/local/prepare_concordance'
include { MERGE_PREPARE_CONCORDANCE             } from '../../modules/local/merge_prepare_concordance'
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
    downsample_rate

    main:
    FASTQ_QC_TRIM_FASTQ_FASTP(input_fastqs)
    FASTQ_ALIGN_BWA(reference, FASTQ_QC_TRIM_FASTQ_FASTP.out.trimmed_reads, bwaidx)
    BAM_DOWNSAMPLE_SAMTOOLS(FASTQ_ALIGN_BWA.out.align.map {[it[0], it[1]]}, downsample_rate)
    BAM_IMPUTE_GLIMPSE2(ref_panel, BAM_DOWNSAMPLE_SAMTOOLS.out.align)
    BCFTOOLS_INDEX(BAM_IMPUTE_GLIMPSE2.out.bcf)
    VARCALL_REF_BCFTOOLS(reference, faidx, FASTQ_ALIGN_BWA.out.align.combine(ref_panel))

    ch_validate = VARCALL_REF_BCFTOOLS.out.bcf
    ch_imputed  = BAM_IMPUTE_GLIMPSE2.out.bcf.join(BCFTOOLS_INDEX.out.csi)
    
    BCF_RENAME_BCFTOOLS(ch_validate.join(ch_imputed))
   
    ch_chrs     = ref_panel.map{it[0]}.collect()
   
    PREPARE_CONCORDANCE(ch_chrs, BCF_RENAME_BCFTOOLS.out.imputed_fixed_bcf)
    
    ch_ref_panel_with_index = ref_panel.map{[it[1], it[2]]}.collect()
    ch_i                    = PREPARE_CONCORDANCE.out.imputed
    ch_v                    = ch_validate.map { [it[0], [it[2], it[3]]] }
        .groupTuple()
        .map { sid, files -> [sid, files.flatten()] }
    ch_i_v                  = ch_i.join(ch_v)
    ch_concordance_txt      = PREPARE_CONCORDANCE.out.concordance_list

    GLIMPSE2_CONCORDANCE(
        ch_ref_panel_with_index, 
        ch_i_v, 
        ch_concordance_txt, 
        '1 5 10 20 50 100 200 500 1000 2000 5000 10000  20000 50000 100000 15000', 
        0, 0
        )

    IMPUTATION_ACCURACY_PLOT(GLIMPSE2_CONCORDANCE.out.rsquare_grp)

    MULTIQC(
        FASTQ_QC_TRIM_FASTQ_FASTP.out.fastp_json                       |
        mix(FASTQ_QC_TRIM_FASTQ_FASTP.out.fastqc)                      |
        mix(FASTQ_QC_TRIM_FASTQ_FASTP.out.fastqc_trimmed)              |
        mix(FASTQ_ALIGN_BWA.out.flagstat.map{it[1]})                   |
        mix(FASTQ_ALIGN_BWA.out.mosdepth.map{it[1]})                   |
        mix(BAM_IMPUTE_GLIMPSE2.out.bcfstats_imputed.map{it[1]})       |
        mix(FASTQ_ALIGN_BWA.out.mosdepth_summary.map{it[1]})           |
        mix(GLIMPSE2_CONCORDANCE.out.errors_grp.map{it[1]})            |
        mix(GLIMPSE2_CONCORDANCE.out.errors_spl.map{it[1]})            |
        collect,
        'summary'
    ) 
}