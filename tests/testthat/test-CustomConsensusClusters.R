# Expected values are derived by hand from the table in helper-synthetic.R.
# The expression of a custom cluster in a sample is the sum of the expression
# of the CTSSs that it contains in that sample, among the CTSSs that pass the
# user-facing filter of CustomConsensusClusters() (signal above threshold in at
# least nrPassThreshold samples).  No other CTSS filter is applied.  The score
# of a cluster is the sum of its normalised expression in all samples.  The
# thresholds below are chosen so that the results do not depend on whether
# the comparison with the threshold is strict or not.
#
# CustomConsensusClusters() renames the clusters after their coordinates, so
# the input clusters are named that way from the start.

cc_high <- "chr1:2000-2060:+"  # Contains 2000, 2005 and 2060.
cc_low  <- "chr1:5000-5010:+"  # Contains 5000 and 5010.

custom_clusters <- function() {
  gr <- GenomicRanges::GRanges(c(cc_high, cc_low))
  names(gr) <- as.character(gr)
  gr
}

test_that("custom clusters made of low CTSSs are retained with the default threshold", {
  # Default threshold = 0, nrPassThreshold = 1: all CTSSs with a signal > 0 in
  # at least one sample are used, that is all of them.
  ce <- CustomConsensusClusters(synthetic_CAGEexp(), custom_clusters())

  # Cluster "low": 5000 (TPM: A = 1) and 5010 (TPM: C = 1 / 2 = 0.5).  Summed
  # over samples, their signal is 1 and 0.5, both <= 1, the threshold of the
  # hidden sum-over-samples filter.  5000 sits exactly at that threshold, which
  # was compared strictly; 5010 is clearly below it.
  expect_true(cc_low %in% cc_ranges(ce))
  # Raw counts: A = 1 (5000); B = 0; C = 1 (5010).
  expect_equal(cc_assay(ce, cc_low, "counts"),     c(A = 1, B = 0, C = 1))
  # TPM:        A = 1; B = 0; C = 0.5.
  expect_equal(cc_assay(ce, cc_low, "normalized"), c(A = 1, B = 0, C = 0.5))
  # Score: 1 + 0 + 0.5 = 1.5.
  expect_equal(cc_score(ce, cc_low), 1.5)

  # Cluster "high": 2000, 2005 and 2060.
  # Raw counts: A = 10 + 3 = 13; B = 4; C = 2.
  expect_equal(cc_assay(ce, cc_high, "counts"),     c(A = 13, B = 4, C = 2))
  # TPM:        A = 13; B = 4 * 2 = 8; C = 2 / 2 = 1.
  expect_equal(cc_assay(ce, cc_high, "normalized"), c(A = 13, B = 8, C = 1))
  # Score: 13 + 8 + 1 = 22.
  expect_equal(cc_score(ce, cc_high), 22)

  # Tags in clusters: A = 13 + 1 = 14; B = 4 + 0 = 4; C = 2 + 1 = 3.
  expect_equal(setNames(ce$outOfClusters, sampleLabels(ce)),
               c(A = 1e6 - 14, B = 5e5 - 4, C = 2e6 - 3))

  # With threshold = 0.75 (TPM), 5000 (A = 1 > 0.75) is used but 5010
  # (C = 0.5) is not.  Raw counts and TPM: A = 1; B = 0; C = 0.  Score: 1.
  ce <- CustomConsensusClusters(synthetic_CAGEexp(), custom_clusters(),
                                threshold = 0.75)
  expect_equal(cc_assay(ce, cc_low, "counts"),     c(A = 1, B = 0, C = 0))
  expect_equal(cc_assay(ce, cc_low, "normalized"), c(A = 1, B = 0, C = 0))
  expect_equal(cc_score(ce, cc_low), 1)
})

test_that("the threshold and nrPassThreshold arguments of CustomConsensusClusters() still apply", {
  # These expectations also pass before the removal of the hidden CTSS filter,
  # because all the CTSSs that they keep have a summed TPM > 1.

  # Among the CTSSs of the custom clusters, those with at least 2 raw tags
  # (threshold = 1.5) in at least 1 sample are 2000 (A = 10), 2005 (A = 3) and
  # 2060 (B = 4).  Not 5000 nor 5010 (1 tag at most), so "low" is empty and
  # is not reported.
  ce <- CustomConsensusClusters(synthetic_CAGEexp(), custom_clusters(),
                                threshold = 1.5, nrPassThreshold = 1,
                                thresholdIsTpm = FALSE)
  expect_identical(cc_ranges(ce), cc_high)
  # Raw counts: A = 10 + 3 = 13; B = 4; C = 2.  TPM: A = 13; B = 8; C = 1.
  expect_equal(cc_assay(ce, cc_high, "counts"),     c(A = 13, B = 4, C = 2))
  expect_equal(cc_assay(ce, cc_high, "normalized"), c(A = 13, B = 8, C = 1))

  # Among the CTSSs of the custom clusters, only 2005 has a TPM above 0.25 in
  # at least 2 samples (A = 3 and C = 1); the others have a signal in a single
  # sample.  Raw counts: A = 3; B = 0; C = 2.  TPM: A = 3; B = 0; C = 1.
  ce <- CustomConsensusClusters(synthetic_CAGEexp(), custom_clusters(),
                                threshold = 0.25, nrPassThreshold = 2)
  expect_identical(cc_ranges(ce), cc_high)
  expect_equal(cc_assay(ce, cc_high, "counts"),     c(A = 3, B = 0, C = 2))
  expect_equal(cc_assay(ce, cc_high, "normalized"), c(A = 3, B = 0, C = 1))
})
