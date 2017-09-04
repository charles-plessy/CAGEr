#' @include AllClasses.R
#' @include CAGEexp.R
#' 
#' @name CAGEr-class
#' 
#' @title CAGEr objects
#' 
#' @details The CAGEr package provides two classes of objects to load, contain
#' and process CAGE data:
#' \itemize{
#'   \item The \code{\link{CAGEset}} class is the original object format in CAGEr,
#'   as when published in Haberle V, Forrest ARR, Hayashizaki Y, Carninci P and
#'   Lenhard B (2015). \dQuote{CAGEr: precise TSS data retrieval and high-resolution
#'   promoterome mining for integrative analyses.} \emph{Nucleic Acids Research},
#'   43, pp. e51., \href{http://nar.oxfordjournals.org/content/43/8/e51}{http://nar.oxfordjournals.org/content/43/8/e51}. 
#'   
#'   \item The \code{\link{CAGEexp}} class is a new class format in 2017, which
#'   is based on the \code{\link{MultiAssayExperiment}} class.  In comparison with
#'   \code{CAGEset}, objects, \code{CAGEexp} objects benefit from a a more efficient
#'   data storage, using \code{DataFrame}s of run-length-encoded (\code{Rle})
#'   integers, allowing for the loading and use of much larger transcriptome datasets.
#' }
#' Most CAGEr functions support both classes interchangabely, and the \code{CAGEr}
#' class was created for methods that support both classes identically.
#' 
#' @docType class
#' @import methods
#' @import BiocGenerics
#' @exportClass CAGEr

setClassUnion("CAGEr", c("CAGEset", "CAGEexp"))


#' @name getRefGenome
#' 
#' @title Attempt to load a BSgenome
#' 
#' @details Internal function that retreives a BSgenome object or throws an error if not available.
#' 
#' @return A BSgenome object
#' 
#' @param reference.genome
#' 
#' @author Charles Plessy
#' 
#' @noRd

getRefGenome <- function(reference.genome) {
  if (is.null(reference.genome))
    stop("Can not run this function with a NULL genome; see ", sQuote('help("genomeName")'), ".")
  if(reference.genome %in% rownames(installed.packages()) == FALSE)
    stop("Requested genome is not installed! Please install required BSgenome package before running CAGEr.")
  requireNamespace(reference.genome)
  getExportedValue(reference.genome, reference.genome)
}

#' @name sampleLabels
#' 
#' @title Extracting CAGE datasets labels from CAGEr objects
#' 
#' @description Extracts the labels and colors of CAGE datasets
#' from \code{\link{CAGEset}} and \code{\link{CAGEexp}} objects.
#' 
#' @param object A CAGEr object.
#' 
#' @return Returns a named character vector of labels of all CAGE datasets
#' present in the CAGEr object.  The values are the lables and the names
#' are the colors.
#' 
#' @note If no colors are supplied, then default colors will be assigned
#' usign the \code{rainbow} function.  Assigned colors are not guaranteed
#' to be stable.
#' 
#' @details Renaming samples is possible only in \code{CAGEexp} objects, before
#' data is loaded.
#' 
#' @author Vanja Haberle
#' 
#' @examples 
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' sampleLabels(exampleCAGEset)
#' 
#' @family CAGEr accessor methods
#' @seealso \code{\link{setColors}}
#' 
#' @importFrom grDevices rainbow
#' @export

setGeneric("sampleLabels", function(object) standardGeneric("sampleLabels"))

#' @rdname sampleLabels

setMethod("sampleLabels", "CAGEset", function (object)
  object@sampleLabels)

#' @rdname sampleLabels

setMethod("sampleLabels", "CAGEexp", function (object){
  sl <- object$sampleLabels
  if (! is.null(object$Colors)) {
    names(sl) <- object$Colors }
  else {
    names(sl) <- rainbow(length(sl))
  }
  sl
})

#' @name validSamples
#' @noRd
#' @title Private function
#' @details Check if a vector of strings or numbers can be used to identify a sample.

setGeneric("validSamples", function(object, x) standardGeneric("validSamples"))

setMethod("validSamples", "CAGEr", function (object, x){
  objName <- deparse(substitute(object))
  if(is.null(x))
      return(TRUE)
  if(class(x) == "character")
    if (all(x %in% sampleLabels(object)))
      return(TRUE)
  if(class(x) == "numeric")
    if (all(x %in% seq_along(sampleLabels(object))))
      return(TRUE)
  stop("Sample(s) not found! Check ", sQuote(paste0("sampleLabels(", objName, ")")), ".")
})
