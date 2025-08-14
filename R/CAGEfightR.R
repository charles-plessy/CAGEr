# Wrappers to CAGEfightR functions.

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
  se <- CTSStagCountSE(object)
  colData(se) <- colData(object)
  rowRanges(se) <- as(rowRanges(se), "StitchedGPos")
  colData(se)$Name <- colData(se)$sampleLabels
  assays(se) <- List(counts=as(as.matrix(as.data.frame(assay(se))), "dgCMatrix"))
  enhancers <- quickEnhancers(se)
  c(enhancers = enhancers, object)
})

# Export normalized CTSS for CAGEfightR
export_normalized_ctss <- function(object = "CAGEexp"){
  # object is copied and modified from CAGEr: the original code only takes the counts
  # Extract CTSS count matrix as SummarizedExperiment
  se <- CTSStagCountSE(object)
  # Clean up and structure metadata
  colData(se) <- colData(object)
  rowRanges(se) <- as(rowRanges(se), "StitchedGPos")
  colData(se)$Name <- colData(se)$sampleLabels

  # Convert counts to sparse matrix to save memory
  assays(se, withDimnames=FALSE) <- List(
      counts = as(as.matrix(as.data.frame(assays(se)[[1]])), "dgCMatrix"),
      TPM = as(as.matrix(as.data.frame(assays(se)[[2]])), "dgCMatrix"))

  # Save as main working object
  cfSampleCTSSs <- se

  return (cfSampleCTSSs)
}
