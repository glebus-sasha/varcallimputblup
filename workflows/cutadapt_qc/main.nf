// Include processes
include { FASTQC as FASTQC_BEFORE } from '../../modules/fastqc'
include { FASTP                   } from '../../modules/fastp'
include { CUTADAPT                } from '../../modules/cutadapt'
include { FASTQC as FASTQC_AFTER  } from '../../modules/fastqc'
include { COVERAGE_SUMMARY_MULTIQC} from '../../modules/cutudapt_summary_multiqc'

workflow CUTADAPT_QC{
    take:
    input_fastqs
    adapters
    
    main:
    FASTQC_BEFORE(input_fastqs)
    CUTADAPT(input_fastqs, adapters)
    FASTQC_AFTER(CUTADAPT.out.cutadapted_reads)
    COVERAGE_SUMMARY_MULTIQC(
        FASTQC_BEFORE.out.zip.collect(),
        FASTQC_AFTER.out.zip.collect(),
}