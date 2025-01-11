// Include processes
include { BAM_BREADTH       } from '../../modules/local/bam_breadth'
include { COV_STATS         } from '../../modules/local/cov_stats'
include { COV_SUMMARY       } from '../../modules/local/cov_summary'
include { GENOME_LENGTH     } from '../../modules/local/genome_length'

workflow COVERAGE_SUMMARY{
    take:
    align
    bcfstats
    mosdepth_summary
    reference

    main:
    align |
    BAM_BREADTH
    GENOME_LENGTH(reference)
    breadth = BAM_BREADTH.out.breadth
    DEPTH_BREADTH(mosdepth_summary.join(breadth), GENOME_LENGTH.out.genome_length)

//    COV_STATS(breadth.join(depth).join(bcfstats))
//    COV_SUMMARY(COV_STATS.out.cov_stats.map{it -> it[1]}.collect())

}