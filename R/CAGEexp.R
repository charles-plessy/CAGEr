#' @rdname CAGEexp
#' 
#' @examples 
#' 
#' pathsToInputFiles <- list.files( system.file("extdata", package = "CAGEr")
#'                                , "ctss$"
#'                                , full.names = TRUE)
#' 
#' myCAGEexp <-
#'   new( "CAGEexp"
#'      , colData = DataFrame( inputFiles = pathsToInputFiles
#'                           , sampleLabels = sub(".chr17.ctss"
#'                                               , ""
#'                                               , basename(pathsToInputFiles)))
#'      , metadata = list( genomeName = "BSgenome.Drerio.UCSC.danRer7"
#'                       , inputFilesType = "ctss"))
#' 
#' getCTSS(myCAGEexp)
#' 
#' colData(myCAGEexp)
#' 
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