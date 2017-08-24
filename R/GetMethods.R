#' @include AggregationFunctions.R ClusteringFunctions.R CAGEr.R

################################################################
# Functions for retrieving data from CAGEset and CAGEexp objects

#' @name genomeName
#' 
#' @title Extracting genome name from CAGEr objects
#' 
#' @description Extracts the name of a referent genome from a
#' \code{\link{CAGEset}} and \code{\link{CAGEexp}} objects.
#' 
#' @param object A CAGEset or CAGEexp object.
#' 
#' @return Returns a name of a BSgenome package used as a referent genome.
#' 
#' @details \code{\link{CAGEexp}} objects constructed with \code{NULL} in place
#' of the genome name can not run some commands that need access to genomic data,
#' such as BigWig export or G-correction.
#' 
#' @family CAGEr accessor methods
#' 
#' @author Vanja Haberle
#' 
#' @examples 
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' genomeName(exampleCAGEset)
#' 
#' @docType methods
#' @rdname genomeName
#' @export

setGeneric(
name="genomeName",
def=function(object){
	standardGeneric("genomeName")
})

setMethod("genomeName",
signature(object = "CAGEset"),
function (object){
	object@genomeName
})

setMethod("genomeName",
signature(object = "CAGEexp"),
function (object){
 metadata(object)$genomeName
})

#' @name inputFiles
#' 
#' @title Extracting paths to input files from CAGEr objects
#' 
#' @description Extracts the paths to CAGE data input files from
#' \code{\link{CAGEset}} and \code{\link{CAGEexp}} objects.
#' 
#' @param object A CAGEset or CAGEexp object.
#' 
#' @return Returns a character vector of paths to CAGE data input files.
#' 
#' @family CAGEr accessor methods
#' 
#' @author Vanja Haberle
#' 
#' @examples 
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' inputFiles(exampleCAGEset)
#' 
#' @docType methods
#' @export

setGeneric(
name="inputFiles",
def=function(object){
	standardGeneric("inputFiles")
})

setMethod("inputFiles",
signature(object = "CAGEset"),
function (object){
	object@inputFiles
})

setMethod("inputFiles",
signature(object = "CAGEexp"),
function (object){
  object$inputFiles
})

#' @name inputFilesType
#'
#' @title Input file formats for CAGEr objects
#' 
#' @description Get or set the information on the type of CAGE data input
#' files from \code{\link{CAGEset}} and \code{\link{CAGEexp}} objects.
#' 
#' @param object A CAGEset or CAGEexp object.
#' 
#' @return Returns the type of the file format of CAGE data input files,
#' \emph{e.g.} \code{"bam"} or \code{"ctss"}.  In the case of \code{CAGEexp}
#' objects, the return value is character vector with one member per sample.
#' 
#' @details The following input file types are supported:
#' 
#' \itemize{
#'   \item{\code{bam}:}
#'   {A single-ended BAM file.}
#' 
#'   \item{\code{bamPairedEnd}:}
#'   {A paired-ended BAM file.}
#' 
#'   \item{\code{bed}:}
#'   {A BED file where each line counts for one molecule.}
#' 
#'   \item{\code{bedScore}:}
#'   {A BED file where the score indicates a number of counts for a
#'   given alignment}
#' 
#'   \item{\code{CAGEscanMolecule}:}
#'   {Experimental.  For the CAGEscan 3.0 pipeline.}
#' 
#'   \item{\code{ctss}:}
#'   {A tabulation-delimited file describing CAGE Transcription
#'   Start Sites (CTSS) with four columns indicating \emph{chromosome},
#'   \emph{1-based coordinate}, \emph{strand} and \emph{score} respectively.}
#' 
#'   \item{\code{CTSStable}}{}
#' 
#'   \item{\code{FANTOM5}}{}
#' 
#'   \item {\code{ENCODE}}{}
#' 
#'   \item{\code{FANTOM3and4}}{}
#' 
#'   \item{\code{ZebrafishDevelopment}}{}
#' }
#' 
#' @author Vanja Haberle
#' 
#' @seealso \code{\link{getCTSS}}
#' 
#' @examples 
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' inputFilesType(exampleCAGEset)
#' 
#' @family CAGEr accessor methods
#' @docType methods
#' @export

setGeneric(
name="inputFilesType",
def=function(object){
	standardGeneric("inputFilesType")
})

setMethod("inputFilesType",
signature(object = "CAGEset"),
function (object){
	object@inputFilesType
})

setMethod("inputFilesType",
signature(object = "CAGEexp"),
function (object){
  object$inputFilesType
})


#' @name librarySizes
#' 
#' @title Extracting library sizes from CAGEr objects
#' 
#' @description Extracts the library sizes (total number of CAGE tags) for all CAGE datasets
#' from \code{\link{CAGEset}} and \code{\link{CAGEexp}} objects.
#' 
#' @param object A CAGEset or CAGEexp object.
#' 
#' @details Library sizes are calculated when loading data with the \code{getCTSS}
#' function and stored in the \code{librarySizes} slot of \code{CAGEset} objects,
#' or in the \code{librarySizes} column of the \code{colData} ov \code{CAGEexp} objects.
#' 
#' @return Returns an integer vector of total number of CAGE tags (library size) for all CAGE
#' datasets in the CAGEr object.
#' 
#' @seealso \code{\link{getCTSS}}
#' 
#' @examples 
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' librarySizes(exampleCAGEset)
#'
#' @author Vanja Haberle
#' @family CAGEr accessor methods
#' @docType methods
#' @export

setGeneric(
name="librarySizes",
def=function(object){
	standardGeneric("librarySizes")
})

setMethod("librarySizes",
signature(object = "CAGEset"),
function (object){
	object@librarySizes
})

setMethod("librarySizes",
signature(object = "CAGEexp"),
function (object){
  as.integer(object$librarySizes)
})

#' @name CTSScoordinates
#' 
#' @title Extracting genomic coordinates of TSSs from CAGEr object
#' 
#' @description Extracts the genomic coordinates of all detected TSSs
#' from \code{\link{CAGEset}} and \code{\link{CAGEexp}} objects.
#' 
#' @param object A CAGEset or CAGEexp object.
#' 
#' @return \code{CTSScoordinates} returns a \code{data.frame} with genomic coordinates of all
#' TSSs. \code{pos} column contains 1-based coordinate of the TSS.
#' 
#' @seealso
#' \code{\link{getCTSS}}
#' 
#' @examples
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' CTSS <- CTSScoordinates(exampleCAGEset)
#' head(CTSS)
#'
#' @author Vanja Haberle
#' @family CAGEr accessor methods
#' @docType methods
#' @export

setGeneric(
name="CTSScoordinates",
def=function(object){
	standardGeneric("CTSScoordinates")
})

setMethod("CTSScoordinates",
signature(object = "CAGEset"),
function (object){
	object@CTSScoordinates
})

setMethod("CTSScoordinates",
signature(object = "CAGEexp"),
function (object){
  gr <- rowRanges(experiments(object)$tagCountMatrix)
  data.frame( chr = as.character(seqnames(gr))
            , pos = start(gr)
            , strand = as.character(strand(gr))
            , stringsAsFactors = FALSE)
})

#' @name CTSScoordinatesGR
#' @rdname CTSScoordinates
#' 
#' @return \code{CTSScoordinatesGR} returns the coordinates as genomic ranges.  A
#' \code{filteredCTSSidx} column metadata will be present if \code{\link{clusterCTSS}}
#' was ran earlier.
#' 
#' @seealso \code{\link{clusterCTSS}}
#' 
#' @author Charles Plessy
#' 
#' @examples
#' CTSScoordinatesGR(exampleCAGEset)
#' 
#' @export

setGeneric(
name="CTSScoordinatesGR",
def=function(object){
	standardGeneric("CTSScoordinatesGR")
})

setMethod("CTSScoordinatesGR",
signature(object = "CAGEset"),
function (object){
  ctssCoord <- object@CTSScoordinates
  ctssCoord <- GRanges(ctssCoord$chr, IRanges(ctssCoord$pos, ctssCoord$pos), ctssCoord$strand)
  genome(ctssCoord) <- object@genomeName
  if(!identical(object@filteredCTSSidx, logical()))
    ctssCoord$filteredCTSSidx <- object@filteredCTSSidx
  ctssCoord
})

setMethod("CTSScoordinatesGR",
signature(object = "CAGEexp"),
function (object){
  rowRanges(CTSStagCountSE(object))
})

#' @name CTSStagCount
#' 
#' @title Extract raw CAGE TSSs expression tables from CAGEr objects
#' 
#' @description Extracts the tag count for all detected TSSs in all CAGE datasets
#'              from \code{\link{CAGEset}} and \code{\link{CAGEexp}} objects.
#' 
#' @param object A CAGEset or CAGEexp object.
#' @param sample (For \code{CTSStagCountGR}.) Name or number of a sample.
#'  
#' @return Returns an object with number of CAGE tags supporting each TSS
#' (rows) in every CAGE dataset (columns).  The class of the object depends on the
#' function being called:
#' 
#' \itemize{
#'   \item{\code{CTSStagCount}:}
#'     {A \code{\link{data.frame}}, containing coordinates and expression values.}
#'   \item{\code{CTSStagCountDf}:}
#'     {A \code{data.frame}, containing only expression values.}
#'   \item{\code{CTSStagCountDF}:}
#'     {A \code{\link{DataFrame}} of \code{\link{Rle}} integers.}
#'   \item{\code{CTSStagCountDA}:}
#'     {A \code{\link{DelayedArray}} wrapping a \code{DataFrame} of \code{Rle} integers.}
#'   \item{\code{CTSStagCountSE}:}
#'     {A \code{\link{RangedSummarizedExperiment}} containing a \code{DataFrame} of
#'      \code{Rle} integers.}
#'   \item{\code{CTSStagCountGR}:}
#'     {A \code{\link{GenomicRanges}} object containing a \code{score} column indicating
#'     expression values for a given sample.}
#' }
#' 
#' @seealso \code{\link{getCTSS}}
#' 
#' @examples
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' head(CTSStagCount(exampleCAGEset))
#' head(CTSStagCountDf(exampleCAGEset))
#' CTSStagCountDF(exampleCAGEset)
#' CTSStagCountDA(exampleCAGEset)
#' 
#' exampleCAGEexp <- readRDS(system.file("extdata", "CAGEexp.rds", package="CAGEr"))
#' head(CTSStagCount(exampleCAGEexp))
#' head(CTSStagCountDf(exampleCAGEset))
#' CTSStagCountDF(exampleCAGEset)
#' CTSStagCountDA(exampleCAGEset)
#' 
#' @author Vanja Haberle
#' @author Charles Plessy
#' 
#' @family CAGEr accessor methods
#' @family CAGEr CTSS methods
#' @docType methods
#' @export

setGeneric(
name="CTSStagCount",
def=function(object){
	standardGeneric("CTSStagCount")
})

setMethod("CTSStagCount",
signature(object = "CAGEset"),
function (object){
	cbind(object@CTSScoordinates, object@tagCountMatrix)
})

setMethod("CTSStagCount",
signature(object = "CAGEexp"),
function (object){
  cbind( CTSScoordinates(object)
       , as.data.frame(lapply(assay(experiments(object)$tagCountMatrix), as.integer)))
})

#' @name CTSStagCountDf
#' @rdname CTSStagCount
#' 
#' @export

setGeneric(
name="CTSStagCountDf",
def=function(object){
	standardGeneric("CTSStagCountDf")
})

setMethod("CTSStagCountDf",
signature(object = "CAGEset"),
function (object){
	object@tagCountMatrix
})

setMethod("CTSStagCountDf",
signature(object = "CAGEexp"),
function (object){
  as.data.frame(lapply(assay(experiments(object)$tagCountMatrix), as.integer))
})

#' @name CTSStagCountDF
#' @rdname CTSStagCount
#'  
#' @export

setGeneric(
name="CTSStagCountDF",
def=function(object){
	standardGeneric("CTSStagCountDF")
})

setMethod("CTSStagCountDF",
signature(object = "CAGEset"),
function (object){
	DF <- object@tagCountMatrix
	DF <- lapply(DF, as.integer)
	DF <- lapply(DF, Rle)
	DataFrame(DF)
})

setMethod("CTSStagCountDF",
signature(object = "CAGEexp"),
function (object){
  assay(CTSStagCountSE(object))
})


#' @name CTSStagCountDA
#' @rdname CTSStagCount
#' 
#' @import DelayedArray DelayedArray
#' @export

setGeneric(
name="CTSStagCountDA",
def=function(object){
	standardGeneric("CTSStagCountDA")
})

setMethod("CTSStagCountDA",
signature(object = "CAGEr"),
function (object){
  DelayedArray(CTSStagCountDF(object))
})

#' @name CTSStagCountGR
#' @rdname CTSStagCount
#'  
#' @export

setGeneric("CTSStagCountGR", function(object, sample) standardGeneric("CTSStagCountGR"))

setMethod( "CTSStagCountGR", "CAGEr", function (object, sample) {
  if (! (sample %in% sampleLabels(object) |
       sample %in% seq_along(sampleLabels(object))))
  stop(sQuote("sample"), " must be the name or number of a sample label.")
  gr <- CTSScoordinatesGR(object)
  score(gr) <- CTSStagCountDF(object)[[sample]]
  gr <- gr[score(gr) != 0]
  gr
})

#' @name CTSStagCountTable
#' 
#' @title Extracting CAGE tag count for TSSs from CAGEr objects
#' 
#' @description Extracts the tag count for all detected TSSs in all CAGE datasets
#' from \code{\link{CAGEset}} and \code{\link{CAGEexp}} objects.
#' 
#' @param object A CAGEset or CAGEexp object.
#'  
#' @return Returns an expression table with the number of CAGE tags supporting each
#' TSS (rows) in every CAGE dataset (columns).  The table is in \code{data.frame}
#' format for \code{CAGEset} objects and in \code{DataFrame} format for \code{CAGEexp}
#' objects.  Use this function when the next consumer can handle both formats.
#' 
#' @seealso \code{\link{getCTSS}}
#' 
#' @examples
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' tagCount <- CTSStagCount(exampleCAGEset)
#' head(tagCount)
#' 
#' @author Vanja Haberle
#' @family CAGEr accessor methods
#' @docType methods
#' @export

setGeneric(
name="CTSStagCountTable",
def=function(object){
	standardGeneric("CTSStagCountTable")
})

setMethod("CTSStagCountTable",
signature(object = "CAGEset"),
  function(object) CTSStagCountDf(object)
)

setMethod("CTSStagCountTable",
signature(object = "CAGEexp"),
  function(object) CTSStagCountDF(object)
)

#' @name CTSStagCountSE
#' @rdname CTSStagCount
#' 
#' @export

setGeneric(
name="CTSStagCountSE",
def=function(object){
	standardGeneric("CTSStagCountSE")
})

setMethod("CTSStagCountSE",
signature(object = "CAGEset"),
function (object){
	  colData <- data.frame(row.names = sampleLabels(object), samplename = sampleLabels(object), samplecolor = names(sampleLabels(object)))
	  if (identical(object@normalizedTpmMatrix, data.frame())) {
	    return(SummarizedExperiment( assays  = list(counts=CTSStagCountDF(object))
	                               , rowData = CTSScoordinatesGR(object)
	                               , colData = colData))
	  }
	  SummarizedExperiment( assays = list( counts              = CTSStagCountDF(object)
	                                     , normalizedTpmMatrix = CTSSnormalizedTpmDF(object))
	                      , rowData = CTSScoordinatesGR(object)
	                      , colData = colData)
})

setMethod("CTSStagCountSE",
signature(object = "CAGEexp"),
function (object){
  experiments(object)$tagCountMatrix
})

#' @name CTSSnormalizedTpm
#' 
#' @title Extracting normalized CAGE signal for TSSs from CAGEr objects
#' 
#' @description Extracts the normalized CAGE signal for all detected TSSs
#' in all CAGE datasets from \code{\link{CAGEset}} and
#' \code{\link{CAGEexp}} objects.
#' 
#' @param object A CAGEset or CAGEexp object.
#' 
#' @return \code{CTSSnormalizedTpm} returns a \code{data.frame} containing coordinates
#' and normalized CAGE signal supporting each TSS (rows) in every CAGE dataset (columns).
#' 
#' @seealso \code{\link{normalizeTagCount}}
#' 
#' @examples 
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' head(CTSSnormalizedTpm(exampleCAGEset))
#' 
#' ce <- readRDS(system.file(package = "CAGEr", "extdata/CAGEexp.rds"))
#' normalizeTagCount(ce)
#' head(CTSSnormalizedTpm(ce))
#' 
#' @author Vanja Haberle
#' @family CAGEr accessor methods
#' @docType methods
#' @export

setGeneric(
name="CTSSnormalizedTpm",
def=function(object){
	standardGeneric("CTSSnormalizedTpm")
})

setMethod("CTSSnormalizedTpm",
signature(object = "CAGEset"),
function (object){
	cbind(object@CTSScoordinates, object@normalizedTpmMatrix)
})

setMethod("CTSSnormalizedTpm",
signature(object = "CAGEexp"),
function (object){
  cbind( CTSScoordinates(object)
       , data.frame(lapply(assays(CTSStagCountSE(ce))[["normalizedTpmMatrix"]], decode)))
})

#' @name CTSSnormalizedTpmDf
#' @rdname CTSSnormalizedTpm
#' 
#' @return \code{CTSSnormalizedTpmDf} returns a \code{data.frame} of normalised expression values.
#' 
#' @examples 
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' head(CTSSnormalizedTpmDf(exampleCAGEset))
#' 
#' ce <- readRDS(system.file(package = "CAGEr", "extdata/CAGEexp.rds"))
#' normalizeTagCount(ce)
#' head(CTSSnormalizedTpmDf(ce))

setGeneric(
name="CTSSnormalizedTpmDf",
def=function(object){
	standardGeneric("CTSSnormalizedTpmDf")
})

setMethod("CTSSnormalizedTpmDf",
signature(object = "CAGEset"),
function (object){
	object@normalizedTpmMatrix
})

setMethod("CTSSnormalizedTpmDf",
signature(object = "CAGEexp"),
function (object){
  data.frame(lapply(assays(CTSStagCountSE(ce))[["normalizedTpmMatrix"]], decode))
})

#' @name CTSSnormalizedTpmDF
#' @rdname CTSSnormalizedTpm
#' 
#' @return \code{CTSSnormalizedTpmDF} returns a \code{DataFrame} of normalised expression values.
#' 
#' @export

setGeneric(
name="CTSSnormalizedTpmDF",
def=function(object){
	standardGeneric("CTSSnormalizedTpmDF")
})

setMethod("CTSSnormalizedTpmDF",
signature(object = "CAGEset"),
function (object){
	DF <- object@normalizedTpmMatrix
	DF <- lapply(DF, Rle)
	DataFrame(DF)
})

setMethod("CTSSnormalizedTpmDF",
signature(object = "CAGEexp"),
function (object){
  assays(object[["tagCountMatrix"]])$normalizedTpmMatrix
})

#' @name CTSSnormalizedTpmGR
#' @rdname CTSSnormalizedTpm
#'  
#' @export

setGeneric("CTSSnormalizedTpmGR", function(object, sample) standardGeneric("CTSSnormalizedTpmGR"))

setMethod( "CTSSnormalizedTpmGR", "CAGEr", function (object, sample) {
  if (! (sample %in% sampleLabels(object) |
       sample %in% seq_along(sampleLabels(object))))
  stop(sQuote("sample"), " must be the name or number of a sample label.")
  gr <- CTSScoordinatesGR(object)
  score(gr) <- CTSSnormalizedTpmDF(object)[[sample]]
  gr <- gr[score(gr) != 0]
  gr
})

#' @name CTSSclusteringMethod
#' 
#' @title Extracting CTSS clustering method from CAGEr objects.
#' 
#' @description Extracts the label of the method used for CTSS clustering into tag
#' clusters from \code{\link{CAGEset}} and \code{\link{CAGEexp}} objects.
#' 
#' @param object A CAGEset or CAGEexp object.
#' 
#' @return Returns a label of the method used for CTSS clustering.
#' 
#' @seealso \code{\link{clusterCTSS}}
#' 
#' @examples 
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' CTSSclusteringMethod(exampleCAGEset)
#' 
#' @author Vanja Haberle
#' @family CAGEr accessor methods
#' @docType methods
#' @export

setGeneric(
name="CTSSclusteringMethod",
def=function(object){
	standardGeneric("CTSSclusteringMethod")
})

setMethod("CTSSclusteringMethod",
signature(object = "CAGEset"),
function (object){
	object@clusteringMethod
})

setMethod("CTSSclusteringMethod",
signature(object = "CAGEexp"),
function (object){
	metadata(object)$clusteringMethod
})


setGeneric("getTagClusterGR", function(object, sample = NULL) {
  validSamples(object, sample)
  standardGeneric("getTagClusterGR")
})

setMethod( "getTagClusterGR", "CAGEset", function (object, sample) {
  if(is.null(sample))
    return(TCdataframe2granges(object@tagClusters))
  TCdataframe2granges(object@tagClusters[[sample]])
})

setMethod( "getTagClusterGR", "CAGEexp", function (object, sample) {
  if(is.null(sample))
    return(metadata(ce)$tagClusters)
  metadata(ce)$tagClusters[[sample]]
})


#' @name tagClusters
#' 
#' @title Extract tag clusters (TCs) for individual CAGE experiments
#' 
#' @description Extracts tag clusters (TCs) produced by \code{\link{clusterCTSS}} function for
#' a specified CAGE experiment from a CAGEr object.
#' 
#' @param object A \code{\link{CAGEr}} object.
#' 
#' @param sample Label of the CAGE dataset (experiment, sample) for which to extract tag clusters.
#' If \code{sample = NULL}, a list of all the clusters for each samples is returned.
#' 
#' @param returnInterquantileWidth Should the interquantile width for each tag cluster be returned.
#' 
#' @param qLow,qUp Position of which quantile should be used as a left (lower) or right (upper)
#' boundary (for \code{qLow} and \code{qUp} respectively) when calculating interquantile width. 
#' Default value \code{NULL} results in using the start coordinate of the cluster.  Used only
#' when \code{returnInterquantileWidth = TRUE}, otherwise ignored.
#' 
#' @return Returns a \code{data.frame} with genomic coordinates, position of dominant TSS, total
#' CAGE signal and additional information for all TCs from specified CAGE dataset (sample).  If
#' \code{returnInterquantileWidth = TRUE}, interquantile width for each TC is also calculated
#' using specified quantile positions and returned in the data frame.
#' 
#' @author Vanja Haberle
#' @family CAGEr accessor methods
#' @family CAGEr clusters functions
#' @export
#' 
#' @examples
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' TC <- tagClusters(exampleCAGEset, "sample2", returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)
#' head(TC)
#' 
#' ce <- readRDS(system.file(package = "CAGEr", "extdata/CAGEexp.rds"))
#' normalizeTagCount(ce)
#' clusterCTSS( object = ce, threshold = 50, thresholdIsTpm = TRUE
#'            , nrPassThreshold = 1, method = "distclu", maxDist = 20
#'            , removeSingletons = TRUE, keepSingletonsAbove = 100)
#' cumulativeCTSSdistribution(ce, clusters = "tagClusters")
#' quantilePositions( ce, clusters = "tagClusters"
#'                  , qLow = c(0.1,0.2), qUp = c(0.8,0.9))
#' TC <- tagClusters(ce, "sample2", TRUE, 0.1, 0.9)
#' head(TC)

setGeneric( "tagClusters"
          , function( object, sample = NULL
                    , returnInterquantileWidth = FALSE, qLow = NULL, qUp = NULL) {
  validSamples(object, sample)
  standardGeneric("tagClusters")
})

setMethod("tagClusters", "CAGEr", function (object, sample, returnInterquantileWidth, qLow, qUp){
  if (is.null(sample)) {
    tc.list <- lapply(sampleLabels(object), tagClusters, object = object, returnInterquantileWidth = returnInterquantileWidth, qLow = qLow, qUp = qUp)
    names(tc.list) <- sampleLabels(object)
    return(tc.list)
  }
  
  if (class(object) == "CAGEset") {
    tc <- object@tagClusters[[sample]]
  } else if (class(object) == "CAGEexp") {
    tc <- TCgranges2dataframe(metadata(object)$tagClusters[[sample]])
  } else {
    stop("Unsupported CAGEr class.")
  }
  
  if (returnInterquantileWidth) {
    if (is.null(qLow) | is.null(qUp))
      stop("No quantiles specified! Please specify which quantile positions should be used to calculate width (qLow and qUp arguments)!")
    if (is.null(tagClustersQuantileLow(object)) & is.null(tagClustersQuantileUp(object)))
      stop("Interquantile width cannot be returned because no quantile positions have been calculated yet! Run 'quantilePositions()' first to get the positions of the desired quantiles!")
    if( (!(paste0("q_", qLow) %in% colnames(tagClustersQuantileLow(object)[[sample]]) &
           paste0("q_", qUp)  %in% colnames(tagClustersQuantileUp(object)[[sample]]))))
      stop("Interquantile width cannot be returned because specified quantile positions have not been calculated! Run 'quantilePositions()' again to get the positions of the desired quantiles!")
  	tc.w <- merge( tagClustersQuantileLow(object)[[sample]]
  	             , tagClustersQuantileUp(object)[[sample]])
  	tc.w <- tc.w[,c("cluster", paste0("q_", qLow), paste0("q_", qUp))]
  	tc.w$interquantile_width <- tc.w[,3] - tc.w[,2] + 1
  	tc <- merge(tc, tc.w)
  }
  tc
})

#' @name filteredCTSSidx
#' @noRd

setGeneric("filteredCTSSidx", function(object) standardGeneric("filteredCTSSidx"))

setMethod("filteredCTSSidx", "CAGEset", function (object){
	object@filteredCTSSidx
})

setMethod("filteredCTSSidx", "CAGEexp", function (object){
  decode(rowData(CTSStagCountSE(ce))$filteredCTSSidx)
})

#' @name tagClustersQuantile
#' @title Quantile metadata stored in CAGEr objects.
#' 
#' @param object A \code{\link{CAGEr}} object.
#' @param value A list (one entry per sample) of data frames with multiple columns:
#'        \code{cluster} for the cluster ID, and then \code{q_0.n} where \code{0.n}
#'        indicates a quantile.
#' @param samples Sample name(s), number(s) or \code{NULL} (default) for all samples.
NULL

#' @name tagClustersQuantileLow
#' @rdname tagClustersQuantile

setGeneric("tagClustersQuantileLow", function(object, samples = NULL) {
  validSamples(object, samples)
  standardGeneric("tagClustersQuantileLow")
})

setMethod("tagClustersQuantileLow", "CAGEset", function (object, samples) {
	if (is.null(samples)) return(object@tagClustersQuantileLow)
  object@tagClustersQuantileLow[[samples]]
})

setMethod("tagClustersQuantileLow", "CAGEexp", function (object, samples) {
  if (is.null(samples)) return(metadata(object)$tagClustersQuantileLow)
  metadata(object)$tagClustersQuantileLow[[samples]]
})


#' @name tagClustersQuantileUp
#' @rdname tagClustersQuantile

setGeneric("tagClustersQuantileUp", function(object, samples = NULL) {
  validSamples(object, samples)
  standardGeneric("tagClustersQuantileUp")
})

setMethod("tagClustersQuantileUp", "CAGEset", function (object, samples) {
  if (is.null(samples)) return(object@tagClustersQuantileUp)
  object@tagClustersQuantileUp[[samples]]
})

setMethod("tagClustersQuantileUp", "CAGEexp", function (object, samples) {
  if (is.null(samples)) return(metadata(object)$tagClustersQuantileUp)
  metadata(object)$tagClustersQuantileUp[[samples]]
})


#' @name getConsensusClusters
#' @noRd

setGeneric("getConsensusClusters", function(object) standardGeneric("getConsensusClusters"))

setMethod("getConsensusClusters", "CAGEset", function (object){
	object@consensusClusters
})

setMethod("getConsensusClusters", "CAGEexp", function (object){
  gr <- rowRanges(consensusClustersSE)
  data.frame()
})

#' @name consensusClustersGR
#' @noRd
#' @export

setGeneric("consensusClustersGR", function(object) standardGeneric("consensusClustersGR"))

setMethod("consensusClustersGR", "CAGEset", function (object)
	CCdataframe2granges(object@consensusClusters))

setMethod("consensusClustersGR", "CAGEexp", function (object)
  rowRanges(consensusClustersSE(object)))

#' @name consensusClustersSE
#' @noRd
#' @export

setGeneric("consensusClustersSE", function(object) standardGeneric("consensusClustersSE"))

setMethod("consensusClustersSE", "CAGEset", function (object)
	stop("CAGEset objects not supported"))

setMethod("consensusClustersSE", "CAGEexp", function (object)
  experiments(object)$consensusClusters)


#' @name consensusClustersDESeq2
#' 
#' @title Export \emph{consensus cluster} expression data for DESeq2 analysis
#' 
#' @description Creates a \code{DESeqDataSet} using the consensus cluster expression
#' data in the experiment slot \code{consensusClusters} and the sample metadata
#' of the \code{\link{CAGEexp}} object.  The formula must be built using factors
#' already present in the sample metadata.
#' 
#' @param object A CAGEexp object.
#' @param design A formula for the DESeq2 analysis.
#' 
#' @author Charles Plessy
#' 
#' @seealso \code{DESeqDataSet} in the \code{DESeq2} package.
#' @family CAGEr expression analysis functions
#' @family CAGEr clusters functions
#' 
#' @examples 
#' ce <- readRDS(system.file(package = "CAGEr", "extdata/CAGEexp.rds"))
#' normalizeTagCount(ce)
#' clusterCTSS( object = ce, threshold = 50, thresholdIsTpm = TRUE
#'            , nrPassThreshold = 1, method = "distclu", maxDist = 20
#'            , removeSingletons = TRUE, keepSingletonsAbove = 100)
#' aggregateTagClusters(ce, tpmThreshold = 50, excludeSignalBelowThreshold = FALSE, maxDist = 100)
#' annotateConsensusClusters(ce, gff)
#' consensusClustersDESeq2(ce, ~group)
#' 
#' @export

setGeneric( "consensusClustersDESeq2"
          , function(object, design) standardGeneric("consensusClustersDESeq2"))

setMethod( "consensusClustersDESeq2", "CAGEset"
         , function (object, design) stop("Not implemented for the CAGEset class."))

setMethod( "consensusClustersDESeq2", "CAGEexp"
         , function (object, design) {
  if (! requireNamespace("DESeq2"))
    stop("This function requires the ", dQuote("DESeq2"), " package; please install it.")
  DESeq2::DESeqDataSetFromMatrix( countData = assay(consensusClustersSE(object))
                                , colData   = colData(object)
                                , rowData   = rowData(consensusClustersSE(object))
                                , design    = design)
})

#' @name consensusClusters
#' 
#' @title Get or set consensus clusters from CAGEr objects
#' 
#' Extracts or inserts the information on consensus clusters from a CAGEr object.
#' 
#' @param object A \code{\link{CAGEr}} object.
#' 
#' @param sample Optional. Label of the CAGE dataset (experiment, sample) for which to extract
#' sample-specific information on consensus clusters.
#' 
#' @param returnInterquantileWidth Should the interquantile width of consensus clusters in
#' specified sample be returned.  Used only when \code{sample} argument is specified, otherwise
#' ignored.
#' 
#' @param qLow Position of which quantile should be used as a left (lower) boundary when
#' calculating interquantile width.  Used only when \code{sample} argument is specified and
#' \code{returnInterquantileWidth = TRUE}, otherwise ignored.
#' 
#' @param qUp Position of which quantile should be used as a right (upper) boundary when
#' calculating interquantile width.  Used only when \code{sample} argument is specified and
#' \code{returnInterquantileWidth = TRUE}, otherwise ignored.
#' 
#' @return Returns a \code{data.frame} with information on consensus clusters, including genomic
#' coordinates.  When \code{sample} argument is NOT specified, total CAGE signal across all CAGE
#' datasets (samples) is returned in the \code{tpm} column.  When \code{sample} argument is
#' specified, the \code{tpm} column contains CAGE signal of consensus clusters in that specific
#' sample.  When \code{returnInterquantileWidth = TRUE}, additional sample-specific information
#' is returned, including position of the dominant TSS, and interquantile width of the consensus
#' clusters in the specified sample.
#' 
#' @author Vanja Haberle
#' 
#' @family CAGEr accessor functions
#' @family CAGEr clusters functions
#' 
#' @examples
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' clusters.general <- consensusClusters(exampleCAGEset)
#' head(clusters.general)
#' 
#' clusters.sample <- consensusClusters(exampleCAGEset, sample = "sample2")
#' head(clusters.sample)
#' 
#' @export

setGeneric(
name="consensusClusters",
def=function(object, sample=NULL, returnInterquantileWidth = FALSE, qLow = NULL, qUp = NULL){
    standardGeneric("consensusClusters")
})

setMethod("consensusClusters",
signature(object = "CAGEr"),
function (object, sample = NULL, returnInterquantileWidth = FALSE, qLow = NULL, qUp = NULL){
	
    if(is.null(sample)){
        return(getConsensusClusters(object))
    }else if(sample %in% sampleLabels(object)){
        
        cc.s <- cbind(cluster = as.integer(rownames(consensusClustersTpm(object))), tpm = consensusClustersTpm(object)[,sample])
        
        if(returnInterquantileWidth & (length(qLow) == 0 | length(qUp) == 0)){
            stop("No quantiles specified! Please specify which quantile positions should be used to calculate width (qLow and qUp arguments)!")
        }else if(returnInterquantileWidth & (length(consensusClustersQuantileLow(object))==0 & length(consensusClustersQuantileUp(object))==0)){
            stop("Interquantile width cannot be returned because no quantile positions for consensus clusters have been calculated yet! Run 'quantilePositions()' first to get the positions of the desired quantiles!")
        }else if(returnInterquantileWidth & (!(paste("q_", qLow, sep = "") %in% colnames(consensusClustersQuantileLow(object)[[sample]]) & paste("q_", qUp, sep = "") %in% colnames(consensusClustersQuantileUp(object)[[sample]])))){
            stop("Interquantile width cannot be returned because specified quantile positions have not been calculated for consensus clusters! Run 'quantilePositions()' again to get the positions of the desired quantiles!")
        }else if(returnInterquantileWidth){
            cc <- getConsensusClusters(object)

            cc.cumsum <- CTSScumulativesConsensusClusters(object)[[sample]]
            a <- lapply(cc.cumsum, function(x) {.get.dominant.ctss(as.numeric(x), isCumulative = T)})
            b <- data.frame(consensus.cluster = as.integer(names(a)), dominant_ctss = unlist(a))
            #cc <- merge(b, cc.s, by.x = 1, by.y = 1, all.x = T, all.y = F)
            
            cc <- merge(cc[,-which(colnames(cc) == "tpm")], b, by.x = 1, by.y = 1, all.x = F, all.y = T)
            cc$dominant_ctss <- cc$start + cc$dominant_ctss
            
            cc <- merge(cc, cc.s, by.x = 1, by.y = 1, all.x = T, all.y = F)
            
            ctss <- CTSSnormalizedTpm(object)[,c("chr", "pos", "strand", sample)]
            cc <- merge(cc, ctss, by.x = c("chr", "strand", "dominant_ctss"), by.y = c("chr", "strand", "pos"), all.x = T, all.y = F)
            colnames(cc)[ncol(cc)] <- "tpm.dominant_ctss"
            cc <- cc[,c("consensus.cluster", "chr", "start", "end", "strand", "dominant_ctss", "tpm", "tpm.dominant_ctss")]
            
            cc.w <- merge(consensusClustersQuantileLow(object)[[sample]], consensusClustersQuantileUp(object)[[sample]])
            cc.w <- cc.w[,c(1, which(colnames(cc.w) == paste("q_", qLow, sep = "")), which(colnames(cc.w) == paste("q_", qUp, sep = "")))]
            cc.w$interquantile_width <- cc.w[,3] - cc.w[,2] + 1
            cc <- merge(cc, cc.w, by.x = "consensus.cluster", by.y = "cluster", all.x = T)
        }else{
            cc <- getConsensusClusters(object)
            cc <- merge(cc[,-which(colnames(cc) == "tpm")], cc.s, by.x = 1, by.y = 1, all.x = T, all.y = F)
            cc <- subset(cc, tpm>0)
        }
        return(cc)
    }else{
        stop("Provided 'sample' not in the CAGE set! Check sampleLabels()")
    }

})


#' @name CTSScumulativesTagClusters
#' 
#' @title Get/set CTSS cumulative TC data
#' 
#' @family CAGEr clusters functions
#' 
#' @export

setGeneric("CTSScumulativesTagClusters", function(object) standardGeneric("CTSScumulativesTagClusters"))

setMethod("CTSScumulativesTagClusters", "CAGEset", function (object){
	object@CTSScumulativesTagClusters
})

setMethod("CTSScumulativesTagClusters", "CAGEexp", function (object){
  metadata(ce)$CTSScumulativesTagClusters
})

#' @name consensusClustersTpm
#' @aliases consensusClustersTpm consensusClustersTpmDf consensusClustersTpmDF
#' 
#' @title Extracting consensus clusters tpm matrix from CAGEset object
#' 
#' @usage
#' consensusClustersTpm(object)
#' consensusClustersTpmDf(object)
#' consensusClustersTpmDF(object)
#' 
#' @description Extracts a table with normalized CAGE tag values for consensus
#' clusters across all samples from a \code{\link{CAGEr}} object.
#' 
#' @param object A CAGEr object.
#' 
#' @return Returns normalized expression valuse of CAGE clusters across all samples.
#' \itemize{
#'   \item For \code{consensusClustersTpm}, a \code{matrix}.
#'   \item For \code{consensusClustersTpmDf}, a \code{data.frame}.
#'   \item For \code{consensusClustersTpmDF}, a \code{DataFrame}.
#' }
#' 
#' @author Vanja Haberle
#' 
#' @family CAGEr clustering methods
#' @seealso \code{\link{consensusClusters}}
#' 
#' @examples
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' clusters.tpm <- consensusClustersTpm(exampleCAGEset)
#' head(clusters.tpm)
#' consensusClustersTpmDF(exampleCAGEset)
#' 
#' @export consensusClustersTpm

setGeneric(
name="consensusClustersTpm",
def=function(object){
	standardGeneric("consensusClustersTpm")
})

setMethod("consensusClustersTpm",
signature(object = "CAGEset"),
function (object){
	object@consensusClustersTpmMatrix
})

setMethod("consensusClustersTpm",
signature(object = "CAGEexp"),
function (object){
  consensusClustersSE <- object[["consensusClusters"]]
  as.matrix(assays(consensusClustersSE)[["normalized"]])
})

#' @export consensusClustersTpmDf

setGeneric(
name="consensusClustersTpmDf",
def=function(object){
	standardGeneric("consensusClustersTpmDf")
})

setMethod("consensusClustersTpmDf",
signature(object = "CAGEset"),
function (object){
	as.data.frame(object@consensusClustersTpmMatrix)
})

setMethod("consensusClustersTpmDf",
signature(object = "CAGEexp"),
function (object){
  consensusClustersSE <- object[["consensusClusters"]]
  as.data.frame(assays(consensusClustersSE)[["normalized"]])
})

#' @export consensusClustersTpmDF

setGeneric(
name="consensusClustersTpmDF",
def=function(object){
	standardGeneric("consensusClustersTpmDF")
})

setMethod("consensusClustersTpmDF",
signature(object = "CAGEset"),
function (object){
	DataFrame(object@consensusClustersTpmMatrix)
})

setMethod("consensusClustersTpmDF",
signature(object = "CAGEexp"),
function (object){
  consensusClustersSE <- object[["consensusClusters"]]
  assays(consensusClustersSE)[["normalized"]]
})

#' expressionClasses
#' @noRd
#' @export

setGeneric(
name="expressionClasses",
def=function(object, what){
	standardGeneric("expressionClasses")
})

setMethod("expressionClasses",
signature(object = "CAGEset"),
function (object, what){
	if(what == "CTSS"){
		classes <- unique(sort(object@CTSSexpressionClasses))
		if(length(classes)>0){
			return(classes)
		}else{
			stop("No expression clustering of CTSSs has been done yet!")
		}
	}else if(what == "consensusClusters") {
		classes <- unique(sort(object@consensusClustersExpressionClasses))
		if(length(classes)>0){
			return(classes)
		}else{
			stop("No expression clustering of consensus clusters has been done yet!")
		}
	}else{
		stop("'what' parameter must be one of the (\"CTSS\", \"consensusClusters\")")
	}
})

#' GeneExpSE
#' 
#' Retreives the SummarizedExperiment containing gene expression levels.
#' 
#' @family CAGEr gene expression analysis functions
#' 
#' @noRd
#' @export

setGeneric(
name="GeneExpSE",
def=function(object){
	standardGeneric("GeneExpSE")
})

setMethod("GeneExpSE",
signature(object = "CAGEset"),
function (object){
	stop("Not implemented for the CAGEset class.")
})

setMethod("GeneExpSE",
signature(object = "CAGEexp"),
function (object){
  experiments(object)$geneExpMatrix
})

#' @name GeneExpDESeq2
#' 
#' @title Export gene expression data for DESeq2 analysis
#' 
#' @description Creates a \code{DESeqDataSet} using the gene expression
#' data in the experiment slot \code{geneExpMatrix} and the sample metadata
#' of the \code{\link{CAGEexp}} object.  The formula must be built using factors
#' already present in the sample metadata.
#' 
#' @param object A CAGEexp object.
#' @param design A formula for the DESeq2 analysis.
#' 
#' @author Charles Plessy
#' 
#' @seealso \code{DESeqDataSet} in the \code{DESeq2} package.
#' @family CAGEr gene expression analysis functions
#' 
#' @examples 
#' ce <- readRDS(system.file(package = "CAGEr", "extdata/CAGEexp.rds"))
#' annotateCTSS(ce, readRDS(system.file("extdata/Zv9_annot.rds", package = "CAGEr")))
#' CTSStoGenes(ce)
#' ce$group <- factor(c("a", "a", "b", "b", "a"))
#' GeneExpDESeq2(ce, ~group)
#' 
#' @export

setGeneric( "GeneExpDESeq2"
          , function(object, design) standardGeneric("GeneExpDESeq2"))

setMethod( "GeneExpDESeq2", "CAGEset"
         , function (object, design) stop("Not implemented for the CAGEset class."))

setMethod( "GeneExpDESeq2", "CAGEexp"
         , function (object, design) {
  if (! requireNamespace("DESeq2"))
    stop("This function requires the ", dQuote("DESeq2"), " package; please install it.")
  DESeq2::DESeqDataSetFromMatrix( countData = assay(GeneExpSE(object))
                                , colData   = colData(object)
                                , rowData   = rowData(GeneExpSE(object))
                                , design    = design)
})


#' turn into tagClustersSE
#' 
#' 
#' #' CTSStagCountSE
#' #' 
#' #' Same as CTSStagCount, but as SummarizedExperiment
#' #' 
#' #' @noRd
#' #' @export
#' 
#' setGeneric(
#' name="CTSStagCountSE",
#' def=function(object){
#' 	standardGeneric("CTSStagCountSE")
#' })
#' 
#' setMethod("CTSStagCountSE",
#' signature(object = "CAGEset"),
#' function (object){
#' 	  colData <- data.frame(row.names = sampleLabels(object), samplename = sampleLabels(object), samplecolor = names(sampleLabels(object)))
#'     SummarizedExperiment(assays = list(counts=CTSStagCountDF(object)), rowData = CTSScoordinatesGR(object), colData = colData)
#' })
#' 
#' setMethod("CTSStagCountSE",
#' signature(object = "CAGEexp"),
#' function (object){
#'   experiments(object)$tagCountMatrix
#' })