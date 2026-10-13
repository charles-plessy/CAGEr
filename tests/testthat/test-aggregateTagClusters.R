# Expected values are derived by hand from the table in helper-synthetic.R and
# from the following definitions.
#
# - Consensus clusters (CCs) are built by merging the tag clusters (TCs) of all
#   samples that have a score >= tpmThreshold, when they overlap or are
#   separated by a gap of less than maxDist bases.
# - The expression of a CC in a sample is the sum of the expression of the
#   CTSSs that it contains in that sample, among the CTSSs that pass the
#   pre-clustering filter (filteredCTSSidx()).  The expression of these CTSSs
#   in other samples, or the number of samples, play no role.
# - The score of a CC is the sum of its normalised expression in all samples.
# - outOfClusters is the library size minus the raw counts in CCs.
#
# All tests use tpmThreshold = 5 and maxDist = 100.  Tag clusters are made by
# distclu(maxDist = 20, keepSingletonsAbove = 0), which joins the CTSSs that
# have a normalised signal > 0 in a sample when they are at most 20 bases
# apart.  Per sample (TPM scores; "pf" = with the pre-clustering filter, which
# removes 1005, 5000 and 5010):
#
#   A: chr1:1000-1015 (2 + 1 + 2 + 2 = 7; pf: 2 + 2 + 2 = 6),
#      chr1:2000-2005 (10 + 3 = 13), chr1:3000-3010 (10 + 3 = 13),
#      chr1:4000-4010 (10 + 3 = 13), chr1:5000 (1; pf: absent),
#      chr2:1000000 (999953)
#   B: chr1:1015 (2), chr1:2060 (8), chr2:1000000 (499995 * 2 = 999990)
#   C: chr1:2005 (1), chr1:2150 (2), chr1:3010 (3), chr1:4010 (6),
#      chr1:5010 (0.5; pf: absent), chr2:1000000 (1999975 / 2 = 999987.5)
#
# TCs with score >= 5:
#   A: 1000-1015, 2000-2005, 3000-3010, 4000-4010, chr2:1000000
#   B: 2060, chr2:1000000
#   C: 4010, chr2:1000000
# Merging them:
#   2000-2005 (A) and 2060 (B): gap of 2059 - 2006 + 1 = 54 bases < 100,
#     hence chr1:2000-2060.  The TC of C at 2150 would have been merged as
#     well (gap of 89 bases), but its score is 2 < 5.
#   4000-4010 (A) and 4010 (C) overlap, hence chr1:4000-4010.
#   1015 (B, score 2 < 5) is not used, but is inside 1000-1015 (A) anyway.
#
# Consensus clusters, with or without the pre-clustering filter:

cc_X <- "chr1:1000-1015:+"  # Broad cluster made only of low CTSSs.
cc_Y <- "chr1:2000-2060:+"
cc_Z <- "chr1:3000-3010:+"
cc_W <- "chr1:4000-4010:+"
cc_P <- "chr2:1000000:+"    # Padding.

aggregate_synthetic <- function(ce, ...)
  aggregateTagClusters(ce, tpmThreshold = 5, maxDist = 100, ...)

test_that("a broad cluster of low CTSSs is retained and quantified from the filteredCTSSidx() CTSSs", {
  # In cluster X, the TC of sample A has a score of 6 (pf) or 7 >= 5, but no
  # CTSS has a normalised signal above 2 in any sample, and the sum over
  # samples is at most 2 + 2 = 4 (position 1015), all <= tpmThreshold = 5.

  # With the pre-clustering filter, 1005 is excluded although it is inside X.
  ce <- aggregate_synthetic(synthetic_tag_clusters(prefilter = TRUE))
  expect_true(cc_X %in% cc_ranges(ce))
  # Raw counts: A = 2 (1000) + 2 (1010) + 2 (1015) = 6; B = 1 (1015); C = 0.
  expect_equal(cc_assay(ce, cc_X, "counts"),     c(A = 6, B = 1, C = 0))
  # TPM:        A = 2 + 2 + 2 = 6; B = 1 * 2 = 2; C = 0.
  expect_equal(cc_assay(ce, cc_X, "normalized"), c(A = 6, B = 2, C = 0))
  # Score: 6 + 2 + 0 = 8.
  expect_equal(cc_score(ce, cc_X), 8)

  # Without the pre-clustering filter, all four CTSSs of X are used.
  ce <- aggregate_synthetic(synthetic_tag_clusters(prefilter = FALSE))
  expect_true(cc_X %in% cc_ranges(ce))
  # Raw counts: A = 2 + 1 + 2 + 2 = 7; B = 1; C = 0.
  expect_equal(cc_assay(ce, cc_X, "counts"),     c(A = 7, B = 1, C = 0))
  # TPM:        A = 7; B = 2; C = 0.
  expect_equal(cc_assay(ce, cc_X, "normalized"), c(A = 7, B = 2, C = 0))
  # Score: 7 + 2 + 0 = 9.
  expect_equal(cc_score(ce, cc_X), 9)
})

test_that("CTSSs that are low in every sample count in their consensus cluster", {
  # Cluster Y contains 2005 (TPM: A = 3, B = 0, C = 1; sum = 4 <= 5).
  ce <- aggregate_synthetic(synthetic_tag_clusters(prefilter = TRUE))
  # Raw counts: A = 10 (2000) + 3 (2005) = 13; B = 4 (2060); C = 2 (2005).
  expect_equal(cc_assay(ce, cc_Y, "counts"),     c(A = 13, B = 4, C = 2))
  # TPM:        A = 10 + 3 = 13; B = 4 * 2 = 8; C = 2 / 2 = 1.
  expect_equal(cc_assay(ce, cc_Y, "normalized"), c(A = 13, B = 8, C = 1))
  # Score: 13 + 8 + 1 = 22.
  expect_equal(cc_score(ce, cc_Y), 22)
})

test_that("expression in a sample does not depend on the number of samples", {
  # Same data, with samples A, B and C, or with samples A and B only.  The
  # consensus clusters are the same, since the TCs of C that pass the
  # threshold (4010 and chr2:1000000) are inside CCs made by A anyway.
  ce3 <- aggregate_synthetic(synthetic_tag_clusters(c("A", "B", "C")))
  ce2 <- aggregate_synthetic(synthetic_tag_clusters(c("A", "B")))
  expect_setequal(cc_ranges(ce3), cc_ranges(ce2))

  # Cluster Z: 3010 has a TPM of 3 in A and 3 in C, never > 5 in one sample;
  # summed over samples it is 3 (A, B) or 3 + 3 = 6 (A, B, C).
  # Raw counts and TPM: A = 10 (3000) + 3 (3010) = 13; B = 0.
  # Cluster W: 4010 has a TPM of 3 in A and 6 in C, > 5 in C only.
  # Raw counts and TPM: A = 10 (4000) + 3 (4010) = 13; B = 0.
  for (ce in list(ce3, ce2)) {
    expect_equal(cc_assay(ce, cc_Z, "counts")[c("A", "B")],     c(A = 13, B = 0))
    expect_equal(cc_assay(ce, cc_Z, "normalized")[c("A", "B")], c(A = 13, B = 0))
    expect_equal(cc_assay(ce, cc_W, "counts")[c("A", "B")],     c(A = 13, B = 0))
    expect_equal(cc_assay(ce, cc_W, "normalized")[c("A", "B")], c(A = 13, B = 0))
  }
  for (assay in c("counts", "normalized"))
    for (cc in cc_ranges(ce2))
      expect_equal(cc_assay(ce3, cc, assay)[c("A", "B")], cc_assay(ce2, cc, assay))

  # Same with excludeSignalBelowThreshold = FALSE (deprecated, see below),
  # which was the only way to see the sum-over-samples filter in isolation.
  ce3 <- suppressWarnings(aggregate_synthetic(synthetic_tag_clusters(c("A", "B", "C")),
                                              excludeSignalBelowThreshold = FALSE))
  ce2 <- suppressWarnings(aggregate_synthetic(synthetic_tag_clusters(c("A", "B")),
                                              excludeSignalBelowThreshold = FALSE))
  expect_equal(cc_assay(ce3, cc_Z, "normalized")[c("A", "B")], c(A = 13, B = 0))
  expect_equal(cc_assay(ce2, cc_Z, "normalized"),              c(A = 13, B = 0))
})

test_that("consensus cluster boundaries are defined by the TCs with score >= tpmThreshold", {
  # Expected from the TC merging derived at the top of this file.  This test
  # also passes before the removal of the hidden CTSS filters, which do not
  # change the boundaries; the presence of cluster X is tested above.
  expected <- c(cc_X, cc_Y, cc_Z, cc_W, cc_P)
  for (prefilter in c(TRUE, FALSE)) {
    ce <- aggregate_synthetic(synthetic_tag_clusters(prefilter = prefilter))
    expect_true(all(cc_ranges(ce) %in% expected))
    expect_true(all(c(cc_Y, cc_Z, cc_W, cc_P) %in% cc_ranges(ce)))
  }
})

test_that("outOfClusters counts the tags outside clusters or outside filteredCTSSidx()", {
  # With the pre-clustering filter, the tags not counted in CCs are:
  #   A: 1005 (1, flagged) + 5000 (1, flagged, outside CCs) = 2
  #   B: none = 0
  #   C: 2150 (4, outside CCs) + 5010 (1, flagged, outside CCs) = 5
  # (Equivalently, for A: 1000000 - (6 + 13 + 13 + 13 + 999953) = 2.)
  ce <- aggregate_synthetic(synthetic_tag_clusters(prefilter = TRUE))
  expect_equal(setNames(ce$outOfClusters, sampleLabels(ce)), c(A = 2, B = 0, C = 5))

  # Without the pre-clustering filter:
  #   A: 5000 (1) = 1;  B: 0;  C: 2150 (4) + 5010 (1) = 5
  ce <- aggregate_synthetic(synthetic_tag_clusters(prefilter = FALSE))
  expect_equal(setNames(ce$outOfClusters, sampleLabels(ce)), c(A = 1, B = 0, C = 5))
})

test_that("excludeSignalBelowThreshold is deprecated and ignored", {
  # A warning naming the argument is emitted when, and only when, the caller
  # supplies it; its wording is not tested.
  ce <- synthetic_tag_clusters()
  expect_no_warning(without <- aggregate_synthetic(ce), message = "excludeSignalBelowThreshold")
  expect_warning(with_true  <- aggregate_synthetic(ce, excludeSignalBelowThreshold = TRUE),
                 "excludeSignalBelowThreshold")
  expect_warning(with_false <- aggregate_synthetic(ce, excludeSignalBelowThreshold = FALSE),
                 "excludeSignalBelowThreshold")
  expect_equal(with_true,  without)
  expect_equal(with_false, without)
})
