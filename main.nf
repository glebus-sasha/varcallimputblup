#!/usr/bin/env nextflow

// Include workflows
include { ALIGN_VARCALL             } from './workflows/align_varcall'
include { IMPUTE                    } from './workflows/impute'

// Include processes
include { MULTIQC                           } from './modules/multiqc'



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

workflow {
    main:
    ALIGN_VARCALL(reference,
        input_fastqs,
        bwaidx,
        faidx)
    IMPUTE(ref_panel_with_index, ref_panel_index, ALIGN_VARCALL.out.align)
    MULTIQC(
        ALIGN_VARCALL.out.fastp.collect(),
        ALIGN_VARCALL.out.fastqc1.collect(),
        ALIGN_VARCALL.out.fastqc2.collect(),
        ALIGN_VARCALL.out.flagstat.collect(),
        ALIGN_VARCALL.out.bcfstats1.collect(),
        IMPUTE.out.bcfstats2.collect()
    )
}



