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
#' @author Charles Plessy
#' @author Katalin Ferenc
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
  # checking content of object on VM
  print("assays counts")
  print(assays(se)$counts)
  print("assays normalizedTpmMatrix")
  print(assays(se)$normalizedTpmMatrix)
  assays(se) <- List(
        counts = as(as.matrix(as.data.frame(assays(se)$counts)), "dgCMatrix"),
        TPM = as(as.matrix(as.data.frame(assays(se)$normalizedTpmMatrix)), "dgCMatrix"))
  # checking content of object on VM
  print(assays(se)$TPM)
  score(rowRanges(se)) <- rowSums(assays(se)$TPM)
  enhancers <- quickEnhancers(se)
  c(enhancers = enhancers, object)
})

