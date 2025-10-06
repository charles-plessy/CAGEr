# Wrappers to CAGEfightR functions.

#' Identify and quantify enhancers.
#' 
#' A convenient wrapper to the function [`CAGEfightR::quickEnhancers()`].
#' 
#' The `CAGEr` object will be converted to a format similar to the output
#' of [`CAGEfightR::quantifyCTSSs()`], and then passed to the `quickEnhancers`
#' function.
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
  colData(se)$Name <- colData(se)$sampleLabels
  assays(se) <- List(
    counts=as(as.matrix(as.data.frame(assays(se)$normalizedTpmMatrix)), "dgCMatrix"))
  score(rowRanges(se)) <- rowSums(assays(se)$TPM)
  enhancers <- quickEnhancers(se)
  c(enhancers = enhancers, object)
})

