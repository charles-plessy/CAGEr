#' CAGEr class to hold all data and metadata about one CAGE experiment.
#' 
#' The \code{CAGEr} class is a \code{\link[MultiAssayExperiment]{MultiAssayExperiment}}
#' object containing all data and metadata about a set of CAGE libraries.  It
#' is a replacement for the \code{\link{CAGEset}} class.  The main difference
#' is that the expression data is stored in \code{\link[S4Vectors]{DataFrame}} objects
#' of \code{\link[S4Vectors]{Rle}}-encoded expression values, instead of plain
#' \code{data.frame}s.  With large datasets, this saves considerable amounts of memory.
#' 
#' @slot metadata A list that must at least contain \code{genomeName} and
#' \code{inputFilesType} members.
#'
#' @import(MultiAssayExperiment)
#' 
#' @examples 
#' 
#' pathsToInputFiles <- list.files( system.file("extdata", package = "CAGEr")
#'                                , "ctss$"
#'                                , full.names = TRUE)
#' sampleLabels <- sub( ".chr17.ctss", "", basename(pathsToInputFiles))
#' myCAGEexp <-
#'   new( "CAGEexp"
#'      , colData = DataFrame( inputFiles = pathsToInputFiles
#'                           , sampleLabels = sampleLabels
#'                           , row.names = sampleLabels)
#'      , metadata = list( genomeName = "BSgenome.Drerio.UCSC.danRer7"
#'                       , inputFilesType = "ctss"))
#' 
#' getCTSS(myCAGEexp)
#' colData(myCAGEexp)
#' 
#' # CTSS data is stored internally as a SummarizedExperiemnt that can be retreived
#' # as a whole, or as GRanges, or as an expression DataFrame.
#' 
#' CTSStagCountSE(myCAGEexp)
#' CTSScoordinatesGR(myCAGEexp)
#' CTSStagCountDF(myCAGEexp)
#' 
#' # Columns of the "colData" table are accessible directly via the "$" operator.
#' 
#' myCAGEexp$librarySizes
#' myCAGEexp$l1 <- colSums(CTSStagCountDf(myCAGEexp) > 0)
#' myCAGEexp$l1
#' 
#' @seealso CAGEset-class
#' 
#' @rdname CAGEexp-class
#' @name CAGEexp-class
#' @export

setClass("CAGEexp",
  contains = "MultiAssayExperiment",
  validity = function(object) {
    #		if(!(object@genomeName %in% suppressWarnings(suppressMessages(BSgenome::available.genomes()))))
    #		return("'genomeName' must be a name of one of the genome packages available in BSgenome! See 'available.genomes()'")
    #		if(object@genomeName %in% rownames(installed.packages()) == FALSE)
    #		return("Requested genome is not installed! Please install required BSgenome package before running CAGEr.")
    
    if (is.null(metadata(object)$genomeName))
      return("Missing BSgenome.")
    
    if (is.null(colData(object)$inputFiles))
      return("Missing input file list.")
    
    if (is.null(metadata(object)$inputFilesType))
      return("Missing input file type.")
    
    if (! metadata(object)$inputFilesType %in%
          c( "bam", "bamPairedEnd", "bed", "ctss", "CTSStable"
           , "FANTOM5", "ENCODE", "FANTOM3and4", "ZebrafishDevelopment"))
      return("'inputFilesType' must be one of supported input file types (\"bam\", \"bamPairedEnd\", \"bed\", \"ctss\", \"CTSStable\")!")
    
    if (is.null(colData(object)$sampleLabels))
      return("Missing sample labels.")
    
    if (!(all(nzchar(colData(object)$sampleLabels))) |
        !(all(substr(colData(object)$sampleLabels, start = 1, stop = 1) %in%c(letters, LETTERS))))
      return("All sample labels must be a non-empty strings beginning with a letter!")
    
    if (length(unique(colData(object)$sampleLabels)) !=
        length(colData(object)$sampleLabels)) 
      return("Duplicated sample labels are not allowed!")
  }
)