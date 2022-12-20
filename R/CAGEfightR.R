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
#' 
#' @examples
#' # The CAGEfightR package needs to be installed before running the command.
#' \dontrun{
#' enhancers <- quickEnhancers(ce)
#' ce <- c(enhancers = enhancers, ce)
#' }
#' 
#' @importFrom CAGEfightR quickEnhancers
#' 
#' @export

setMethod("quickEnhancers", signature(object = "CAGEexp"), function(object) {
  se <- CTSStagCountSE(object)
  colData(se) <- colData(object)
  rowRanges(se) <- as(rowRanges(se), "StitchedGPos")
  colData(se)$Name <- colData(se)$sampleLabels
  assays(se) <- List(counts=as(as.matrix(as.data.frame(assay(se))), "dgCMatrix"))
  quickEnhancers(se)
})
