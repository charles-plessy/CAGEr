# Wrappers to CAGEfightR functions.

#' Helper function to export tag count data into summarized experiment object with StitchedGPos.
#' 
#' @param object A `CAGEexp` object
#' 
#' @return A `SummarizedExperiment` object.
#' 
#' @family CAGEr accessor methods
#' @family CAGEr normalized data functions

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

#' Export of normalized CTSS in an object that can be used as input to CAGEfightR.
#' 
#' An export function for integration of normalized values into CAGEfightR.
#' 
#' @note At the moment the conversion is expensive as it goes from `DataFrame`
#' of `Rle` to `data.frame` to `matrix`.
#' 
#' @param object A `CAGEexp` object with CTSS values
#' 
#' @return A `SummarizedExperiment` object with StitchedGPos values.
#' 
#' @family CAGEr accessor methods
#' @family CAGEr normalized data functions
#' 
#' @examples
#' exportNormalizedCTSS(exampleCAGEexp)
setGeneric("exportNormalizedCTSS", function(object)
  standardGeneric("exportNormalizedCTSS"))


#' @export
#' @rdname exportNormalizedCTSS
#' @aliases exportNormalizedCTSS,CAGEexp-method
#' 
setMethod("exportNormalizedCTSS", signature( object = "CAGEexp"), function(object) {
  se <- .export_tag_counts(object)

  # Convert counts to sparse matrix to save memory
  assays(se, withDimnames=FALSE) <- List(
      counts = as(as.matrix(as.data.frame(assays(se)[[1]])), "dgCMatrix"),
      TPM = as(as.matrix(as.data.frame(assays(se)[[2]])), "dgCMatrix"))

  # Return summarized experiment object
  se
})
