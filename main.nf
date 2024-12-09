#!/usr/bin/env nextflow

// Include processes
include { FASTQC    } from './modules/fastqc'
include { FASTP     } from './modules/fastp'

// Logging pipeline information
log.info """\
\033[0;36m  ==========================================  \033[0m
\033[0;34m       v a r c a l l i m p u t b l u p        \033[0m
\033[0;36m  ==========================================  \033[0m
    """
    .stripIndent(true)

// Define the input channel for reference file
reference = params.reference ? Channel.fromPath("${params.reference}").collect(): null

// Define the input channel for FASTQ files, if provided
input_fastqs = Channel.fromFilePairs(["${params.reads}/*[rR]{1,2}*.*{fastq,fq}*", "${params.reads}/*_{1,2}.{fastq,fq}*"])

// Define the input channel for bwa index files, if provided
bwaidx = params.bwaidx ? Channel.fromPath("${params.bwaidx}/*", checkIfExists: true).collect() : null

// Define the input channel for fai index files, if provided
faidx = params.bwaidx ? Channel.fromPath("${params.faidx}/*.fai", checkIfExists: true).collect() : null

// Define the workflow
workflow test { 
    input_fastqs.view()
}

workflow {
    test()
}



