// Include processes
include { BAM_BREADTH   } from '../../modules/local/breadth'
include { BAM_DEPTH     } from '../../modules/local/depth'
include { COV_STATS     } from '../../modules/local/cov_stats'
include { COV_SUMMARY   } from '../../modules/local/cov_summary'

workflow COVERAGE_SUMMARY{
    take:
    align
    bcfstats

    main:
    align |
    BAM_BREADTH & BAM_DEPTH

    breadth = BAM_BREADTH.out.breadth
    depth = BAM_DEPTH.out.depth_stats

    COV_STATS(breadth.join(depth).join(bcfstats))
    COV_SUMMARY(COV_STATS.out.cov_stats.map{it -> it[1]}.collect())
}