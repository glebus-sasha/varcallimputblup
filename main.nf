#!/usr/bin/env nextflow

// Include processes
include { FASTQC as FASTQC1           } from './modules/fastqc'
include { FASTP                       } from './modules/fastp'
include { FASTQC as FASTQC2           } from './modules/fastqc'
include { BWA_MEM                     } from './modules/bwa_mem'
include { SAMTOOLS_FLAGSTAT           } from './modules/samtools/flagstat'
include { SAMTOOLS_INDEX              } from './modules/samtools/index'
include { BCFTOOLS_MPILEUP            } from './modules/bcftools/mpileup'
include { BCFTOOLS_STATS              } from './modules/bcftools/stats'
include { MULTIQC                     } from './modules/multiqc'
include { GLIMPSE2_CHUNK              } from './modules/glimpse2/chunk'
include { GLIMPSE2_CONCORDANCE        } from './modules/glimpse2/concordance'
include { GLIMPSE2_LIGATE             } from './modules/glimpse2/ligate'
include { GLIMPSE2_PHASE              } from './modules/glimpse2/phase'
include { GLIMPSE2_SPLITREFERENCE     } from './modules/glimpse2/splitreference'

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

ref_panel = Channel.fromPath("${params.ref_panel}").map{file->[file.baseName, file]}
ref_panel_index = Channel.fromPath("${params.ref_panel_index}").map{file->[file.baseName, file]}

bam = Channel.fromPath("${params.bam}/*.bam").map{file->[file.name, file]}
bamindex = Channel.fromPath("${params.bam}/*.bam.bai").map{file->[file.baseName, file]}

// Define the workflow
workflow test { 
    take:
    ref_panel
    ref_panel_index
    bam
    bamindex

    main:
    ref_panel.view()
    ref_panel_index.view()
    GLIMPSE2_CHUNK(ref_panel.join(ref_panel_index))
    IRG_ORG = GLIMPSE2_CHUNK.out.chunk_chr.splitCsv(header:false,sep:'\t').map{T->[T[2],T[3]]}
    GLIMPSE2_SPLITREFERENCE(ref_panel, ref_panel_index, GLIMPSE2_CHUNK.out.chunk_chr, IRG_ORG)
    GLIMPSE2_PHASE(GLIMPSE2_SPLITREFERENCE.out.bin_ref, ref_panel_index, bam.join(bamindex))
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
    BCFTOOLS_STATS(BCFTOOLS_MPILEUP.out.vcf)
    MULTIQC(
        FASTP.out.json.collect(), 
        FASTQC1.out.zip.collect(), 
        FASTQC2.out.zip.collect(), 
        SAMTOOLS_FLAGSTAT.out.flagstat.collect(), 
        BCFTOOLS_STATS.out.bcfstats.collect()
        )
}

workflow IMPUTE {
    take:
    ref_panel

    main:
    GLIMPSE2_CHUNK(ref_panel, ref_panel_index, 'chr1')

    GLIMPSE2_CONCORDANCE
    GLIMPSE2_LIGATE
    GLIMPSE2_PHASE()
    GLIMPSE2_SPLITREFERENCE
}

workflow {
    test(ref_panel, ref_panel_index, bam, bamindex)
}



