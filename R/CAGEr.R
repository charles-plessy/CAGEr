#' @include CAGEexp.R CTSS.R
NULL
 
#' CAGEr objects
#' 
#' The _CAGEr_ package provides one class of objects to load, contain and
#' process CAGE data: the [`CAGEexp`] class, introduced 2017, which is based on the
#' [`MultiAssayExperiment`] class.  In comparison with the original `CAGEset`
#' class (removed in 2021) `CAGEexp` objects benefit from a a more efficient data storage, using
#' `DataFrame`s of run-length-encoded (`Rle`) integers, allowing for the
#' loading and use of much larger transcriptome datasets.
#' 
#' @references Haberle V, Forrest ARR, Hayashizaki Y, Carninci P and Lenhard B
#' (**2015**). \dQuote{CAGEr: precise TSS data retrieval and high-resolution
#' promoterome mining for integrative analyses.} _Nucleic Acids Research_,
#' 43, pp. e51., <http://nar.oxfordjournals.org/content/43/8/e51>
#' 
#' @import methods
#' @import BiocGenerics
#' @exportClass CAGEr

setClassUnion("CAGEr", c("CAGEexp"))  # Legacy of dual support for CAGEset before 2021


#' @name getRefGenome
#' 
#' @title Return a BSgenome or throws an error
#' 
#' @details If the `reference.genome` object exists and is a BSgenome_, it will
#' be returned.  This allows the user to run things like
#' `GenomeInfoDb::seqlevelsStyle(BSgenome.Hsapiens.UCSC.hg19) <- "NCBI"` on
#' the BSgenome object before running `getCTSS`.  If the `reference.genome`
#' object does not exist, attempts to load it and return it, or throws an
#' error if not available.
#' 
#' @return A BSgenome object of the same name as the `reference.genome` argument.
#' 
#' @param reference.genome
#' 
#' @author Charles Plessy
#' @importFrom utils installed.packages
#' 
#' @noRd

getRefGenome <- function(reference.genome) {
  if (is.null(reference.genome))
    stop("Can not run this function with a NULL genome; see ", sQuote('help("genomeName")'), ".")
  if (exists(reference.genome))
    if ("BSgenome" %in% class(get(reference.genome)))
      return(get(reference.genome))
    else
      stop("The ", sQuote(reference.genome), " object in the namespace is not a BSgenome")
  if(reference.genome %in% rownames(installed.packages()) == FALSE)
    stop("Requested genome is not installed! Please install required BSgenome package before running CAGEr.")
  requireNamespace(reference.genome)
  getExportedValue(reference.genome, reference.genome)
}

#' @name sampleLabels
#' 
#' @title Get and set sample labels
#' 
#' @description \code{sampleLabels} gets or sets the labels and colors of CAGE datasets
#' (samples) from \code{\link{CAGEr}} objects.
#' 
#' @param object A CAGEr object.
#' 
#' @return \code{sampleLabels} returns a named character vector representing labels of all
#' CAGE datasets present in the CAGEr object.  The vector values are the labels and the
#' vector names are the colors.
#' 
#' @note If no colors are supplied, then default colors will be assigned
#' usign the \code{rainbow} function.  Assigned colors are not guaranteed
#' to be stable.
#' 
#' @details In \code{CAGEexp} objects, renaming samples is possible only before
#' data is loaded.
#' 
#' @author Vanja Haberle
#' 
#' @examples 
#' sampleLabels(exampleCAGEexp)
#' 
#' @family CAGEr accessor methods
#' @seealso \code{\link{setColors}}
#' 
#' @importFrom grDevices rainbow
#' @export

setGeneric("sampleLabels", function(object) standardGeneric("sampleLabels"))

#' @rdname sampleLabels

setMethod("sampleLabels", "CAGEexp", function (object) {
  sl <- object$sampleLabels
  if (! is.null(object$Colors)) {
    names(sl) <- object$Colors }
  else {
    names(sl) <- rainbow(length(sl))
  }
  sl
})

#' @rdname sampleLabels

setMethod("sampleLabels", "CTSS", function (object)
  object@metadata$sampleLabels)

#' @description \code{sampleList} is an accessory function for convenience
#' iteration in functions such as \code{\link{lapply}} or \code{\link{mapply}}.
#' There is no set method for \code{sampleList}.
#' 
#' @return \code{sampleList} returns a named list where elements and their
#' names are the sample names, for instance: \code{list(sampleA = "sampleA",
#' sampleB = "sampleB")}. Thus, after iterating on it with \code{lapply}, the
#' element names will be sample names.
#' 
#' @examples 
#' sampleList(exampleCAGEexp)
#' 
#' @export
#' @rdname sampleLabels

setGeneric("sampleList", function(object) standardGeneric("sampleList"))

#' @rdname sampleLabels

setMethod("sampleList", "CAGEr", function (object) {
  l <- sampleLabels(object)
  names(l) <- l
  l
})

#' @name validSamples
#' @noRd
#' @title Private function
#' @examples
#' CAGEr:::validSamples(exampleCAGEexp, 1)
#' CAGEr:::validSamples(exampleCAGEexp, "Zf.high")
#' CAGEr:::validSamples(exampleCAGEexp, sampleLabels(exampleCAGEexp))
#' CAGEr:::validSamples(exampleCAGEexp, "all")
#' CAGEr:::validSamples(exampleCAGEexp, NULL)
#' @details Check if a vector of strings or numbers can be used to identify a sample.

setGeneric("validSamples", function(object, x) standardGeneric("validSamples"))

setMethod("validSamples", "CAGEr", function (object, x){
  objName <- deparse(substitute(object))
  if(is.null(x))
      return(TRUE)
  if(inherits(x, "character"))
    if (x == "all")
      return(TRUE)
    if (all(x %in% sampleLabels(object)))
      return(TRUE)
  if(inherits(x, "integer") | inherits(x, "numeric"))
    if (all(x %in% seq_along(sampleLabels(object))))
      return(TRUE)
  stop("Sample(s) not found! Check ", sQuote(paste0("sampleLabels(", objName, ")")), ".")
})


#' Flag CTSSes based on sample expression
#' 
#' Flag CTSSes for that do not pass an expression threshold in at least a given
#' number of samples.  This is typically used to ignore CTSSes that have been
#' seen only once in a single sample, as they can be considered to not be
#' reproduced.
#' 
#' @param object An object from the _CAGEr_ package that contains expression
#' values for multiple samples.
#'
#' @param threshold Flag CTSSs with signal `< threshold`.
#'      
#' @param nrPassThreshold Only flag CTSSs when signal is below threshold in at
#'        least `nrPassThreshold` samples.
#'        
#' @param thresholdIsTpm Logical, is threshold raw tag count value (`FALSE`) or
#'        normalized signal (`TRUE`).
#'        
#' @returns `flagLowExpCTSS` returns a [`Rle`] vector where `TRUE` indicates the
#' index of a CTSS that passes the filter.
#' 
#' @export
#' 
#' @family CAGEr filter functions
#' 
#' @examples
#' flagLowExpCTSS(exampleCAGEexp, threshold = 100, nrPassThreshold = 2)

setGeneric("flagLowExpCTSS",  function( object
                                      , threshold       = 1
                                      , nrPassThreshold = 1
                                      , thresholdIsTpm  = TRUE)
  standardGeneric("flagLowExpCTSS")
)

#' @rdname flagLowExpCTSS

setMethod("flagLowExpCTSS", "CAGEr", function (object, threshold, nrPassThreshold, thresholdIsTpm) {
  flagLowExpCTSS(CTSStagCountSE(object), threshold, nrPassThreshold, thresholdIsTpm)
})

#' @rdname flagLowExpCTSS

setMethod("flagLowExpCTSS", "RangedSummarizedExperiment", function (object, threshold, nrPassThreshold, thresholdIsTpm) {
  assay <- ifelse(thresholdIsTpm, "normalizedTpmMatrix", "counts")
  if(assay == "normalizedTpmMatrix" & is.null(assays(object)[[assay]]))
    stop("Normalise the CAGEr object first with ", sQuote("normalizeTagCount()"), ".")
  flagLowExpCTSS(assays(object)[[assay]], threshold, nrPassThreshold, thresholdIsTpm)
})

#' @rdname flagLowExpCTSS

setMethod("flagLowExpCTSS", "DataFrame", function (object, threshold, nrPassThreshold, thresholdIsTpm) {
  nr.pass.threshold <- rowSums.RleDataFrame(lapply(object, \(x) x > threshold) |> DataFrame())
  nr.pass.threshold >= min(nrPassThreshold, ncol(object))
})

#' @rdname flagLowExpCTSS

setMethod("flagLowExpCTSS", "matrix", function (object, threshold, nrPassThreshold, thresholdIsTpm) {
  nr.pass.threshold <- rowSums(object > threshold)
  nr.pass.threshold >= min(nrPassThreshold, ncol(object))
})

#' @rdname flagLowExpCTSS
#' 
#' @return `filterLowExpCTSS` returns the `CAGEr` object where the output of
#' `flagLowExpCTSS` was stored internally.
#' 
#' @export

setGeneric("filterLowExpCTSS",  function( object
                                        , threshold       = 1
                                        , nrPassThreshold = 1
                                        , thresholdIsTpm  = TRUE)
  standardGeneric("filterLowExpCTSS")
)

#' @rdname flagLowExpCTSS

setMethod("filterLowExpCTSS", "CAGEr", function (object, threshold, nrPassThreshold, thresholdIsTpm) {
  filteredCTSSidx(object) <- flagLowExpCTSS(CTSStagCountSE(object), threshold, nrPassThreshold, thresholdIsTpm)
  object
})
