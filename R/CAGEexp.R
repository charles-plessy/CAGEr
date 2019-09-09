#' CAGEr class to hold all data and metadata about one CAGE experiment.
#' 
#' The [`CAGEr`] class is a [`MultiAssayExperiment`] object containing all data
#' and metadata about a set of CAGE libraries.  It is a replacement for the
#' [`CAGEset`] class.  The main difference is that the expression data is stored
#' in [`DataFrame`] objects of [`Rle`]-encoded expression values, instead of plain
#' `data.frame`s.  With large datasets, this saves considerable amounts of memory.
#' 
#' @slot metadata A list that must at least contain a `genomeName` member.
#' 
#' @details If `genomeName` is `NULL`, checks of chromosome names will be
#' disabled and G-correction will not be possible.  See
#' <https://support.bioconductor.org/p/86437/> for an example on how to create a
#' _BSgenome_ package.
#' 
#' Sample labels must be _syntactically valid_ in the sense of the [make.names()]
#' function, because they will be used as column names in some tables.
#' 
#' @examples 
#' pathsToInputFiles <- list.files( system.file("extdata", package = "CAGEr")
#'                                , "ctss$"
#'                                , full.names = TRUE)
#' sampleLabels <- sub( ".chr17.ctss", "", basename(pathsToInputFiles))
#' 
#' # The CAGEexp object can be created using specific constructor commands
#'                               
#' exampleCAGEexp <-
#'   CAGEexp( genomeName     = "BSgenome.Drerio.UCSC.danRer7"
#'          , inputFiles     = pathsToInputFiles
#'          , inputFilesType = "ctss"
#'          , sampleLabels   = sub( ".chr17.ctss", "", basename(pathsToInputFiles)))
#'          
#' # Alternatively, it can be created just like another MultiAssayExperiment.
#' # This is useful when providing pre-existing colData with many columns.
#' 
#' exampleCAGEexp <-
#'   CAGEexp( metadata = list(genomeName = "BSgenome.Drerio.UCSC.danRer7")
#'          , colData  = DataFrame( inputFiles     = pathsToInputFiles
#'                                , sampleLabels   = sampleLabels
#'                                , inputFilesType = "ctss"
#'                                , row.names      = sampleLabels))
#' 
#' 
#' # Expression data is loaded by the getCTSS() function, that also calculates
#' # library sizes and store them in the object's column data.
#' 
#' getCTSS(exampleCAGEexp)
#' librarySizes(exampleCAGEexp)
#' colData(exampleCAGEexp)
#' 
#' # CTSS data is stored internally as a SummarizedExperiemnt that can be retreived
#' # as a whole, or as GRanges, or as an expression DataFrame.
#' 
#' CTSStagCountSE(exampleCAGEexp)
#' CTSScoordinatesGR(exampleCAGEexp)
#' CTSStagCountDF(exampleCAGEexp)
#' 
#' # Columns of the "colData" table are accessible directly via the "$" operator.
#' 
#' exampleCAGEexp$l1 <- colSums(CTSStagCountDf(exampleCAGEexp) > 0)
#' exampleCAGEexp$l1
#' 
#' @seealso [`CAGEset-class`], [`make.names`]
#' 
#' @rdname CAGEexp-class
#' @aliases CAGEexp-class
#' @import MultiAssayExperiment
#' @import SummarizedExperiment
#' @importFrom BSgenome installed.genomes
#' @importFrom S4Vectors SimpleList
#' @export CAGEexp
#' @exportClass CAGEexp

CAGEexp <- setClass("CAGEexp",
  contains = "MultiAssayExperiment",
  validity = function(object) {
    #		if(!(object@genomeName %in% suppressWarnings(suppressMessages(BSgenome::installed.genomes()()))))
    #		return("'genomeName' must be a name of one of the genome packages available in BSgenome! See 'BSgenome::installed.genomes()'")
    #		if(object@genomeName %in% rownames(installed.packages()) == FALSE)
    #		return("Requested genome is not installed! Please install required BSgenome package before running CAGEr.")
    
    if (! is.null(genomeName(object)))
      if (! genomeName(object) %in% installed.genomes())
        return( paste0(sQuote("genomeName"), " must be the name an installed genome package. "
              , "See ", sQuote("BSgenome::installed.genomes()")
              , " and ", sQuote("BSgenome::available.genomes()"), "."))
    
    if (is.null(colData(object)$inputFiles))
      return("Missing input file list.")
    
    if (is.null(object$inputFilesType))
      return("Missing input file type.")
    
    supportedTypes <- c("bam", "bamPairedEnd", "bed", "bedScore", "bedctss", "CAGEscanMolecule", "ctss", "CTSStable")
    
    if (! all(inputFilesType(object) %in% supportedTypes))
      return( paste(sQuote("inputFilesType"), "must be one of supported input file types:"
            , paste(sQuote(supportedTypes), collapse = ", "), "."))
 
    if (is.null(rownames(colData(object))))
      return("Rownames are missing in colData().  Did you forget to specify them?")
       
    if (is.null(colData(object)$sampleLabels))
      return("Missing sample labels.")
    
    if( !(identical(colData(object)$sampleLabels, make.names(colData(object)$sampleLabels))))
      stop( "Names of samples must be syntactically valid (see "
          , dQuote("?make.names"), ") because they will be used as column names "
          , "in some tables.  In brief, avoid minus signs and do not start with "
          , "numbers.")
    
    if (!(all(nzchar(colData(object)$sampleLabels))) |
        !(all(substr(colData(object)$sampleLabels, start = 1, stop = 1) %in%c(letters, LETTERS))))
      return("All sample labels must be a non-empty strings beginning with a letter!")
    
    if (length(unique(colData(object)$sampleLabels)) !=
        length(colData(object)$sampleLabels)) 
      return("Duplicated sample labels are not allowed!")
  }
)

setMethod( "initialize", "CAGEexp"
         , function( .Object
                   , ...
                   , metadata       = list()
                   , colData        = S4Vectors::DataFrame()
                   , genomeName     = NULL
                   , inputFiles     = NULL
                   , inputFilesType = NULL
                   , sampleLabels   = NULL) {
  if (missing(colData)) {
    if (missing(inputFiles)) stop ("Please provide input files paths.")
    colData <- S4Vectors::DataFrame(inputFiles = inputFiles)
  }
  if (! missing(inputFiles))     colData$inputFiles     <- inputFiles
  if (! missing(inputFilesType)) colData$inputFilesType <- inputFilesType
  if (! missing(sampleLabels))   colData$sampleLabels   <- sampleLabels
  if (! missing(sampleLabels))   rownames(colData)      <- sampleLabels
  if (! missing(genomeName))     metadata$genomeName    <- genomeName
  callNextMethod( .Object
                , metadata = metadata
                , colData  = colData
                , ...)
})

#' Example CAGEexp object.
#' 
#' Lazy-loaded example CAGEexp object, containing most of the CAGEr data
#' structures created with the CAGEr modifier functions.
#' 
#' @format A [CAGEexp] object.
#' 
#' @examples 
#' \dontrun{
#' pathsToInputFiles <- list.files( system.file("extdata", package = "CAGEr")
#'                                , "ctss$"
#'                                , full.names = TRUE)
#' sampleLabels <- sub( ".chr17.ctss", "", basename(pathsToInputFiles))
#' exampleCAGEexp <-
#'   CAGEexp( genomeName     = "BSgenome.Drerio.UCSC.danRer7"
#'          , inputFiles     = pathsToInputFiles
#'          , inputFilesType = "ctss"
#'          , sampleLabels   = sub( ".chr17.ctss", "", basename(pathsToInputFiles)))
#' getCTSS(exampleCAGEexp)
#' librarySizes(exampleCAGEexp)
#' colData(exampleCAGEexp)
#' exampleCAGEexp$l1 <- NULL
#' exampleCAGEexp <- exampleCAGEexp[,c(5, 2, 1, 3, 4)] # Non-aplhabetic order may help catch bugs
#' CTSStagCountSE(exampleCAGEexp) <- CTSStagCountSE(exampleCAGEexp)[1:5000,]  # Slim the object
#' exampleCAGEexp$librarySizes <- sapply(CTSStagCountDF(exampleCAGEexp), sum) # Repair metadata
#' summariseChrExpr(exampleCAGEexp)
#' annotateCTSS(exampleCAGEexp, exampleZv9_annot)
#' CTSStoGenes(exampleCAGEexp)
#' normalizeTagCount(exampleCAGEexp)
#' clusterCTSS(exampleCAGEexp)
#' cumulativeCTSSdistribution(exampleCAGEexp, "tagClusters")
#' quantilePositions(exampleCAGEexp, "tagClusters")
#' aggregateTagClusters(exampleCAGEexp)
#' annotateConsensusClusters(exampleCAGEexp, exampleZv9_annot)
#' cumulativeCTSSdistribution(exampleCAGEexp, "consensusClusters")
#' quantilePositions(exampleCAGEexp, "consensusClusters")
#' save(exampleCAGEexp, file = "data/exampleCAGEexp.RData", compress = "xz")
#' }

"exampleCAGEexp"


#' @name coerce,data.frame,CAGEexp-method
#' @rdname CAGEexp-class

setAs("data.frame", "CAGEexp", function(from){
  
  if((ncol(from) < 4) | !(all(colnames(from)[1:3] == c("chr", "pos", "strand"))))
    stop( "First three columns of the input data.frame must contain chromosome name, "
        , "genomic position and strand of individual TSSs, and must be named 'chr', 'pos' "
        , "and 'strand', respectively!")
  
  if(ncol(from) <4)
    stop( "Input data.frame needs to contain at least one column with CAGE tag counts, "
        , "in addition to first three columns specifying chromosome name, genomic position "
        , "and strand of individual TSSs!")
  
  if(!(is.integer(from[,"pos"])))
    stop( "The 'pos' column in the input data.frame can contain only non-zero integers "
        , "as these are interpreted as 1-based genomic coordinates of TSSs! Make sure the "
        , "'pos' column is of class 'integer'!")
  
  if(any(from[,"pos"] <= 0))
    stop( "The 'pos' column in the input data.frame can contain only non-zero integers as "
        , "these are interpreted as 1-based genomic coordinates of TSSs!")
  
  if(!(all(from[,"strand"] %in% c("+", "-"))))
    stop( "The 'strand' column in the input data.frame can contain only '+' or '-'!")
  
  if(!(all(apply(from[,4:ncol(from),drop=FALSE], 2, is.integer))))
    stop( "The columns specifying CAGE tag counts must be non-negative integers! Make sure "
        , "these columns are of class 'integer'!")

  if(any(apply(from[,4:ncol(from),drop=FALSE], 2, function(x) {any(x < 0)})))
    stop("The columns specifying CAGE tag counts must be non-negative integers!")
  
  sample.labels <- colnames(from)[4:ncol(from)]
  gr            <- GRanges(from$chr, IRanges(from$pos, width = 1), from$strand)
  counts        <- DataFrame(lapply(from[,4:ncol(from),drop=FALSE], Rle))
  
  object <-
    CAGEexp( metadata = list(genomeName = NULL)
           , colData = DataFrame( inputFiles     = "data.frame"
                                , sampleLabels   = sample.labels
                                , inputFilesType = "CTSStable"
                                , row.names      = sample.labels))

  object$librarySizes <- as.integer(colSums(from[,4:ncol(from),drop=FALSE]))
  
  CTSStagCountSE(object) <-
    SummarizedExperiment( rowRanges = gr
                        , assays    = SimpleList(counts = counts))
  
  return(object)
})