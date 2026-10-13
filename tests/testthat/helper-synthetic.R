# Small synthetic CAGEexp objects for testing how consensus clusters are built
# and quantified.
#
# All CTSSs are on the plus strand of chr1, except one "padding" CTSS on chr2
# whose only purpose is to bring each library to a round size, so that the
# "simpleTpm" normalisation gives values that are easy to compute by hand:
#
#   library size   A = 1,000,000  ->  TPM = count
#                  B =   500,000  ->  TPM = count * 2
#                  C = 2,000,000  ->  TPM = count / 2
#
#             Raw counts            Normalised (TPM)
#   pos       A    B    C           A    B    C        Comment
#   1000      2    0    0           2    0    0        Region X: broad cluster
#   1005      1    0    0           1    0    0          made of low CTSSs;
#   1010      2    0    0           2    0    0          1005 is removed by the
#   1015      2    1    0           2    2    0          pre-clustering filter.
#   2000     10    0    0          10    0    0        Region Y: TCs of A and B
#   2005      3    0    2           3    0    1          are merged; the low TC
#   2060      0    4    0           0    8    0          of C at 2150 is not.
#   2150      0    0    4           0    0    2
#   3000     10    0    0          10    0    0        Region Z: 3010 is low in
#   3010      3    0    6           3    0    3          every sample.
#   4000     10    0    0          10    0    0        Region W: 4010 is high in
#   4010      3    0   12           3    0    6          sample C only.
#   5000      1    0    0           1    0    0        Region L: low CTSSs only.
#   5010      0    0    1           0    0    0.5
#   chr2:1000000 (padding)
#        999953 499995 1999975
#
#   Sums without padding: A = 7 + 13 + 13 + 13 + 1  = 47
#                         B = 1 +  4                =  5
#                         C = 2 +  4 +  6 + 12 + 1  = 25
#   Library sizes:        A =   999953 + 47 = 1000000
#                         B =   499995 +  5 =  500000
#                         C =  1999975 + 25 = 2000000

synthetic_ctss_table <- function() {
  data.frame(
    seqnames = c(rep("chr1", 14), "chr2"),
    pos = c(1000L, 1005L, 1010L, 1015L,
            2000L, 2005L, 2060L, 2150L,
            3000L, 3010L,
            4000L, 4010L,
            5000L, 5010L,
            1000000L),
    A   = c(2L, 1L, 2L, 2L,   10L, 3L, 0L, 0L,   10L, 3L,   10L,  3L,   1L, 0L,
            999953L),
    B   = c(0L, 0L, 0L, 1L,    0L, 0L, 4L, 0L,    0L, 0L,    0L,  0L,   0L, 0L,
            499995L),
    C   = c(0L, 0L, 0L, 0L,    0L, 2L, 0L, 4L,    0L, 6L,    0L, 12L,   0L, 1L,
            1999975L)
  )
}

# Build a normalised CAGEexp object from the table above, keeping only the
# samples given as argument.  Like getCTSS(), CTSSs that have no tag in any of
# the kept samples are not included in the object.
#
# filterLowExpCTSS() is run with a threshold that flags no CTSS, because the
# consensusClustersGR() accessor, which aggregateTagClusters() and
# CustomConsensusClusters() call internally, fails when the filteredCTSSidx
# column is absent.  Every CTSS has at least 1 tag, so none is flagged whether
# the comparison with the threshold of 0 is strict or not.

synthetic_CAGEexp <- function(samples = c("A", "B", "C")) {
  tab    <- synthetic_ctss_table()
  counts <- tab[, samples, drop = FALSE]
  keep   <- rowSums(counts) > 0
  tab    <- tab[keep, ]
  counts <- counts[keep, , drop = FALSE]
  stopifnot(all(colSums(counts) %in% c(5e5, 1e6, 2e6)))

  ce <- CAGEexp(
    metadata = list(genomeName = NULL),
    colData  = S4Vectors::DataFrame(
      inputFiles     = paste0(samples, ".ctss"),
      inputFilesType = "ctss",
      sampleLabels   = samples,
      row.names      = samples))

  CTSStagCountSE(ce) <- SummarizedExperiment::SummarizedExperiment(
    rowRanges = CTSS(tab$seqnames, tab$pos, "+"),
    assays    = S4Vectors::SimpleList(
      counts = S4Vectors::DataFrame(lapply(counts, S4Vectors::Rle))))
  ce$librarySizes <- as.integer(colSums(counts))

  ce <- suppressMessages(normalizeTagCount(ce, method = "simpleTpm"))
  filterLowExpCTSS(ce, threshold = 0, nrPassThreshold = 1, thresholdIsTpm = FALSE)
}

# Tag clustering with distclu().  When `prefilter = TRUE`, CTSSs that do not
# have at least 2 raw tags in at least one sample are flagged by
# filterLowExpCTSS() before clustering, so that filteredCTSSidx() is not all
# TRUE.  In the table above, these are 1005, 5000 and 5010.  Otherwise,
# filteredCTSSidx() is TRUE for all CTSSs.  The threshold of 1.5 does not fall
# on a tag count, so the result does not depend on whether the comparison
# with the threshold is strict or not.

synthetic_tag_clusters <- function(samples = c("A", "B", "C"), prefilter = FALSE) {
  ce <- synthetic_CAGEexp(samples)
  if (prefilter)
    ce <- filterLowExpCTSS(ce, threshold = 1.5, nrPassThreshold = 1,
                           thresholdIsTpm = FALSE)
  suppressMessages(distclu(ce, maxDist = 20, keepSingletonsAbove = 0))
}

# Per-sample expression of a consensus cluster, as a named numeric vector.
# Returns NAs if the cluster is missing, so that tests fail on expectations
# instead of stopping with a subscript error.

cc_assay <- function(ce, cluster, assay = c("counts", "normalized")) {
  assay <- match.arg(assay)
  m <- SummarizedExperiment::assay(consensusClustersSE(ce), assay)
  if (! cluster %in% rownames(m))
    return(stats::setNames(rep(NA_real_, ncol(m)), colnames(m)))
  m[cluster, , drop = TRUE]
}

# Score of a consensus cluster as stored in the consensus clusters experiment
# (the consensusClustersGR() accessor recomputes it instead).  NA if the
# cluster is missing.

cc_score <- function(ce, cluster) {
  gr <- SummarizedExperiment::rowRanges(consensusClustersSE(ce))
  if (! cluster %in% names(gr)) return(NA_real_)
  as.numeric(gr[cluster]$score)
}

# Consensus cluster coordinates as character strings, for example
# "chr1:2000-2060:+".

cc_ranges <- function(ce) {
  unname(as.character(SummarizedExperiment::rowRanges(consensusClustersSE(ce))))
}
