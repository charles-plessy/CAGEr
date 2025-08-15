# Wrappers to CAGEfightR functions.

# Helper function to export tag count data
.export_tag_counts <- function(object = "CAGEexp") {
  se <- CTSStagCountSE(object)
  colData(se) <- colData(object)
  rowRanges(se) <- as(rowRanges(se), "StitchedGPos")
  colData(se)$Name <- colData(se)$sampleLabels
  se
}

#' Identify and quantify enhancers.
#' 
#' A convenient wrapper to the function [`CAGEfightR::quickEnhancers()`].
#' 
#' The `CAGEr` object will be converted to a format similar to the output
#' of [`CAGEfightR::quantifyCTSSs()`], and then passed to the `quickEnhancers`
#' function.
#' 
#' @note At the moment the conversion is expensive as it goes from `DataFrame`
#' of `Rle` to `data.frame` to `matrix`.
#' 
#' @param object A `CAGEexp` object
#' 
#' @return A `RangedSummarizedExperiment` object.  See the example below on
#' how to attach it to the experiment list of a `CAGEexp` object.
#' 
#' @family CAGEfightR
#' @family CAGEr object modifiers
#' 
#' @examples
#' # Can not run as long as the test data has nothing on the minus strand!
#' \dontrun{
#' quickEnhancers(exampleCAGEexp)
#' }
#' 
#' @importFrom CAGEfightR quickEnhancers

setGeneric("quickEnhancers", function(object)
  standardGeneric("quickEnhancers"))

#' @export
#' @rdname quickEnhancers
#' @aliases quickEnhancers,CAGEexp-method

setMethod("quickEnhancers", signature(object = "CAGEexp"), function(object) {
  se <- .export_tag_counts(object)
  assays(se) <- List(counts=as(as.matrix(as.data.frame(assay(se))), "dgCMatrix"))
  enhancers <- quickEnhancers(se)
  c(enhancers = enhancers, object)
})


#' Enhancer calling
#'
#' @param ce CAGEexp object with CTSS values
#' @param cfBalanceThreshold threshold for the cagefightr balance score
#' @param unexpressed threshold above which normalized CTSS are considered expressed
#' @param minSamples non inlcusive lower threshold for number of samples supporting enhancers (i.e. where there is bidirectionality)
#' @return enhancers
#' @examples
#' cagefightr_enhancers(
#' ce,
#' cfBalanceThreshold = 0.95,
#' unexpressed = 0,
#' minSamples = 0
#' )
setMethod("CAGEfightREnhancers"
        , signature( object = "CAGEexp", cfBalanceThreshold
                    , unexpressed, minSamples)
                    , function(object) {
  se <- .export_tag_counts(object)

  # Convert counts to sparse matrix to save memory
  assays(se, withDimnames=FALSE) <- List(
      counts = as(as.matrix(as.data.frame(assays(se)[[1]])), "dgCMatrix"),
      TPM = as(as.matrix(as.data.frame(assays(se)[[2]])), "dgCMatrix"))

  # Save as main working object
  cfSampleCTSSs <- se

  # Calculate pooled signal across all samples (average TPM)
  cfSampleCTSSs <- CAGEfightR::calcPooled(
      cfSampleCTSSs,
      inputAssay = "TPM")

  # Calculate how many samples support expression at each CTSS
  cfSampleCTSSs <- CAGEfightR::calcSupport(
      cfSampleCTSSs,
      inputAssay = "counts",
      outputColumn = "support",
      unexpressed = 0)

  # Find bidirectional clusters (potential enhancers)
  sampleBCs <- CAGEfightR::clusterBidirectionally(
      cfSampleCTSSs,
      balanceThreshold = cfBalanceThreshold)

  # Filter bidirectional clusters that are supported in at least 1 sample
  finalSampleBCs <- CAGEfightR::subsetByBidirectionality(
      sampleBCs,
      samples = cfSampleCTSSs,
      minSamples = 0)

  finalSampleBCs
})
