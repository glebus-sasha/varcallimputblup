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

// Define the input channel for reference file
reference = Channel.fromPath("${params.reference}").collect()

// Define the input channel for FASTQ files, if provided
input_fastqs = Channel.fromFilePairs(["${params.reads}/*[rR]{1,2}*.*{fastq,fq}*", "${params.reads}/*_{1,2}.{fastq,fq}*"], flat = true)

// Define the input channel for bwa index files, if provided
bwaidx = Channel.fromPath("${params.bwaidx}/*", checkIfExists: true).collect()

// Define the input channel for fai index files, if provided
faidx = Channel.fromPath("${params.faidx}/*.fai", checkIfExists: true).collect()

// Define the workflow
workflow test { 
    input_fastqs.view()
}

workflow FASTQ_QC_TRIM_ALIGN_VARCALL { 
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
        SAMTOOLS_FLAGSTAT.out.flagstat.collect(), 
        BCFTOOLS_STATS.out.bcfstats.collect()
        )
}

workflow IMPUTE {
    GLIMPSE2_CHUNK
    GLIMPSE2_CONCORDANCE
    GLIMPSE2_LIGATE
    GLIMPSE2_PHASE
    GLIMPSE2_SPLITREFERENCE
}

workflow {
    FASTQ_QC_TRIM_ALIGN_VARCALL()
}



