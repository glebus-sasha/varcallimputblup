#!/usr/bin/env nextflow

// Include processes
include { FASTQC as FASTQC1                 } from './modules/fastqc'
include { FASTP                             } from './modules/fastp'
include { FASTQC as FASTQC2                 } from './modules/fastqc'
include { BWA_MEM                           } from './modules/bwa_mem'
include { SAMTOOLS_FLAGSTAT                 } from './modules/samtools/flagstat'
include { SAMTOOLS_INDEX                    } from './modules/samtools/index'
include { BCFTOOLS_MPILEUP                  } from './modules/bcftools/mpileup'
include { BCFTOOLS_INDEX                    } from './modules/bcftools/index'
include { BCFTOOLS_STATS as BCFTOOLS_STATS1 } from './modules/bcftools/stats'
include { MULTIQC                           } from './modules/multiqc'
include { GLIMPSE2_CHUNK                    } from './modules/glimpse2/chunk'
include { GLIMPSE2_CONCORDANCE              } from './modules/glimpse2/concordance'
include { GLIMPSE2_LIGATE                   } from './modules/glimpse2/ligate'
include { GLIMPSE2_PHASE                    } from './modules/glimpse2/phase'
include { GLIMPSE2_SPLITREFERENCE           } from './modules/glimpse2/splitreference'
include { BCFTOOLS_STATS as BCFTOOLS_STATS2 } from './modules/bcftools/stats'

// Logging pipeline information
log.info """\
\033[0;36m  ==========================================  \033[0m
\033[0;34m       v a r c a l l i m p u t b l u p        \033[0m
\033[0;36m  ==========================================  \033[0m
    """
    .stripIndent(true)

reference = Channel.fromPath("${params.reference}").collect()
input_fastqs = Channel.fromFilePairs(["${params.reads}/*[rR]{1,2}*.*{fastq,fq}*", "${params.reads}/*_{1,2}.{fastq,fq}*"], flat: true)
bwaidx = Channel.fromPath("${params.bwaidx}/*", checkIfExists: true).collect()
faidx = Channel.fromPath("${params.faidx}/*.fai", checkIfExists: true).collect()

ref_panel = Channel.fromPath("${params.ref_panel}").map{file->[file.simpleName, file]}
ref_panel_index = Channel.fromPath("${params.ref_panel_index}").map{file->[file.simpleName, file]}
ref_panel_with_index = ref_panel.join(ref_panel_index)

bam = Channel.fromPath("${params.bam}/*.bam").map{file->[file.simpleName, file]}
bamindex = Channel.fromPath("${params.bam}/*.bam.bai").map{file->[file.simpleName, file]}
align = bam.join(bamindex)

workflow test { 

}

workflow FASTQ_QC_TRIM_ALIGN_VARCALL { 
    take:
    reference
    input_fastqs
    bwaidx
    faidx

    main:
    FASTQC1(input_fastqs)
    FASTP(input_fastqs)
    FASTQC2(FASTP.out.trimmed_reads)
    BWA_MEM(FASTP.out.trimmed_reads, reference, bwaidx)
    SAMTOOLS_FLAGSTAT(BWA_MEM.out.bam)
    SAMTOOLS_INDEX(BWA_MEM.out.bam)
    BCFTOOLS_MPILEUP(reference, SAMTOOLS_INDEX.out.bai, faidx)
    BCFTOOLS_INDEX(BCFTOOLS_MPILEUP.out.bcf)
    BCFTOOLS_STATS1(BCFTOOLS_MPILEUP.out.bcf)
    MULTIQC(
        FASTP.out.json.collect(),
        FASTQC1.out.zip.collect(),
        FASTQC2.out.zip.collect(),
        SAMTOOLS_FLAGSTAT.out.flagstat.collect(),
        BCFTOOLS_STATS1.out.bcfstats.collect()
        )
    emit:
    align = SAMTOOLS_INDEX.out.bai
    bcfstats = BCFTOOLS_STATS1.out.bcfstats
}

workflow BCF_IMPUTE {
    take:
    ref_panel_with_index
    ref_panel_index
    align

    main:
    GLIMPSE2_CHUNK(ref_panel_with_index)
    IRG_ORG = GLIMPSE2_CHUNK.out.chunk_chr.splitCsv(header:false,sep:'\t').map{coord->[coord[2],coord[3]]}
    GLIMPSE2_SPLITREFERENCE(ref_panel_with_index.combine(IRG_ORG))
    GLIMPSE2_PHASE(
        align.combine(GLIMPSE2_SPLITREFERENCE.out.bin_ref.map{it->it[1]}).combine(ref_panel_index)
        )
    GLIMPSE2_LIGATE(GLIMPSE2_PHASE.out.phased_variants.groupTuple())
    BCFTOOLS_STATS2(GLIMPSE2_LIGATE.out.merged_variants)
}

workflow {
    FASTQ_QC_TRIM_ALIGN_VARCALL(reference,
        input_fastqs,
        bwaidx,
        faidx)
    BCF_IMPUTE(ref_panel_with_index, ref_panel_index, FASTQ_QC_TRIM_ALIGN_VARCALL.out.align)
}



