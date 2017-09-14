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
#' genomeName(exampleCAGEset)
#' 
#' @export

setGeneric(
name="genomeName",
def=function(object){
	standardGeneric("genomeName")
})

#' @rdname genomeName

setMethod("genomeName",
signature(object = "CAGEset"),
function (object){
	object@genomeName
})

#' @rdname genomeName

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
#' inputFiles(exampleCAGEset)
#' 
#' @export

setGeneric(
name="inputFiles",
def=function(object){
	standardGeneric("inputFiles")
})

#' @rdname inputFiles

setMethod("inputFiles",
signature(object = "CAGEset"),
function (object){
	object@inputFiles
})

#' @rdname inputFiles

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
#' inputFilesType(exampleCAGEset)
#' 
#' @family CAGEr accessor methods
#' @export

setGeneric(
name="inputFilesType",
def=function(object){
	standardGeneric("inputFilesType")
})

#' @rdname inputFilesType

setMethod("inputFilesType",
signature(object = "CAGEset"),
function (object){
	object@inputFilesType
})

#' @rdname inputFilesType

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
#' librarySizes(exampleCAGEset)
#'
#' @author Vanja Haberle
#' @family CAGEr accessor methods
#' @export

setGeneric(
name="librarySizes",
def=function(object){
	standardGeneric("librarySizes")
})

#' @rdname  librarySizes

setMethod("librarySizes",
signature(object = "CAGEset"),
function (object){
	object@librarySizes
})

#' @rdname  librarySizes

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
#' CTSS <- CTSScoordinates(exampleCAGEset)
#' head(CTSS)
#'
#' @author Vanja Haberle
#' @family CAGEr accessor methods
#' @export

setGeneric(
name="CTSScoordinates",
def=function(object){
	standardGeneric("CTSScoordinates")
})

#' @rdname CTSScoordinates

setMethod("CTSScoordinates",
signature(object = "CAGEset"),
function (object){
	object@CTSScoordinates
})

#' @rdname CTSScoordinates

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
#' @importFrom GenomeInfoDb genome genome<-
#' @importFrom IRanges IRanges
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

#' @rdname CTSScoordinates

setMethod("CTSScoordinatesGR",
signature(object = "CAGEset"),
function (object){
  ctssCoord <- object@CTSScoordinates
  ctssCoord <- GRanges(ctssCoord$chr, IRanges(ctssCoord$pos, ctssCoord$pos), ctssCoord$strand)
  genome(ctssCoord) <- object@genomeName
  if(!identical(object@filteredCTSSidx, logical()))
    ctssCoord$filteredCTSSidx <- object@filteredCTSSidx
  if(!identical(CTSSexpressionClasses(object), character()))
    if (length(CTSSexpressionClasses(object)) == length(ctssCoord)) {
      ctssCoord$CTSSexpressionClasses <- CTSSexpressionClasses(object)
    } else {
      warning("Skipping expression classes: not same length as CTSS.")
    }
  ctssCoord
})

#' @rdname CTSScoordinates

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
#' @param samples (For \code{CTSStagCountGR}.) Name(s) or number(s) identifying sample(s).
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
#'     {A \code{CTSS} object (wrapping \code{GRanges}) containing a \code{score}
#'     column indicating expression values for a given sample.}
#' }
#' 
#' @seealso \code{\link{getCTSS}}
#' 
#' @examples
#' head(CTSStagCount(exampleCAGEset))
#' head(CTSStagCountDf(exampleCAGEset))
#' CTSStagCountDF(exampleCAGEset)
#' CTSStagCountDA(exampleCAGEset)
#' 
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
#' @export

setGeneric(
name="CTSStagCount",
def=function(object){
	standardGeneric("CTSStagCount")
})

#' @rdname CTSStagCount

setMethod("CTSStagCount",
signature(object = "CAGEset"),
function (object){
	cbind(object@CTSScoordinates, object@tagCountMatrix)
})

#' @rdname CTSStagCount

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

#' @rdname CTSStagCount

setMethod("CTSStagCountDf",
signature(object = "CAGEset"),
function (object){
	object@tagCountMatrix
})

#' @rdname CTSStagCount

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

#' @rdname CTSStagCount

setMethod("CTSStagCountDF",
signature(object = "CAGEset"),
function (object){
	DF <- object@tagCountMatrix
	DF <- lapply(DF, as.integer)
	DF <- lapply(DF, Rle)
	DataFrame(DF)
})

#' @rdname CTSStagCount

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

#' @rdname CTSStagCount

setMethod("CTSStagCountDA",
signature(object = "CAGEr"),
function (object){
  DelayedArray(CTSStagCountDF(object))
})

#' @rdname CTSStagCount


#' @name CTSStagCountGR
#' @rdname CTSStagCount
#'  
#' @export

setGeneric("CTSStagCountGR", function(object, samples) standardGeneric("CTSStagCountGR"))

#' @rdname CTSStagCount

setMethod( "CTSStagCountGR", "CAGEr", function (object, samples) {
  if (! (samples %in% sampleLabels(object) |
       samples %in% seq_along(sampleLabels(object))))
  stop(sQuote("samples"), " must be the name or number of a sample label.")
  gr <- CTSScoordinatesGR(object)
  score(gr) <- CTSStagCountDF(object)[[samples]]
  gr <- gr[score(gr) != 0]
  .CTSS(gr)
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
#' tagCount <- CTSStagCount(exampleCAGEset)
#' head(tagCount)
#' 
#' @author Vanja Haberle
#' @family CAGEr accessor methods
#' @export

setGeneric(
name="CTSStagCountTable",
def=function(object){
	standardGeneric("CTSStagCountTable")
})

#' @rdname CTSStagCountTable

setMethod("CTSStagCountTable",
signature(object = "CAGEset"),
  function(object) CTSStagCountDf(object)
)

#' @rdname CTSStagCountTable

setMethod("CTSStagCountTable",
signature(object = "CAGEexp"),
  function(object) CTSStagCountDF(object)
)


#' @name CTSStagCountSE
#' @rdname CTSStagCount
#' 
#' @export

setGeneric("CTSStagCountSE", function(object) standardGeneric("CTSStagCountSE"))

#' @rdname CTSStagCountTable

setMethod("CTSStagCountSE", "CAGEset", function (object) {
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

#' @rdname CTSStagCountTable

setMethod("CTSStagCountSE", "CAGEexp", function (object) {
  se <- experiments(object)$tagCountMatrix
  if (is.null(se)) stop("Could not find CTSS tag counts, see ", sQuote("?getCTSS"), ".")
  se
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
#' head(CTSSnormalizedTpm(exampleCAGEset))
#' 
#' normalizeTagCount(exampleCAGEexp)
#' head(CTSSnormalizedTpm(exampleCAGEexp))
#' 
#' @author Vanja Haberle
#' @family CAGEr accessor methods
#' @export

setGeneric(
name="CTSSnormalizedTpm",
def=function(object){
	standardGeneric("CTSSnormalizedTpm")
})

#' @rdname CTSSnormalizedTpm

setMethod("CTSSnormalizedTpm",
signature(object = "CAGEset"),
function (object){
	cbind(object@CTSScoordinates, object@normalizedTpmMatrix)
})

#' @rdname CTSSnormalizedTpm

setMethod("CTSSnormalizedTpm",
signature(object = "CAGEexp"),
function (object){
  cbind( CTSScoordinates(object)
       , data.frame(lapply(assays(CTSStagCountSE(object))[["normalizedTpmMatrix"]], decode)))
})


#' @name CTSSnormalizedTpmDf
#' @rdname CTSSnormalizedTpm
#' 
#' @return \code{CTSSnormalizedTpmDf} returns a \code{data.frame} of normalised expression values.
#' 
#' @examples 
#' head(CTSSnormalizedTpmDf(exampleCAGEset))
#' 
#' normalizeTagCount(exampleCAGEexp)
#' head(CTSSnormalizedTpmDf(exampleCAGEexp))
#' 
#' @export

setGeneric(
name="CTSSnormalizedTpmDf",
def=function(object){
	standardGeneric("CTSSnormalizedTpmDf")
})

#' @rdname CTSSnormalizedTpm

setMethod("CTSSnormalizedTpmDf",
signature(object = "CAGEset"),
function (object){
	object@normalizedTpmMatrix
})

#' @rdname CTSSnormalizedTpm

setMethod("CTSSnormalizedTpmDf",
signature(object = "CAGEexp"),
function (object){
  data.frame(lapply(assays(CTSStagCountSE(object))[["normalizedTpmMatrix"]], decode))
})


#' @name CTSSnormalizedTpmDF
#' @rdname CTSSnormalizedTpm
#' 
#' @return \code{CTSSnormalizedTpmDF} returns a \code{DataFrame} of normalised expression values.
#' 
#' @importFrom SummarizedExperiment assays
#' 
#' @export

setGeneric(
name="CTSSnormalizedTpmDF",
def=function(object){
	standardGeneric("CTSSnormalizedTpmDF")
})

#' @rdname CTSSnormalizedTpm

setMethod("CTSSnormalizedTpmDF",
signature(object = "CAGEset"),
function (object){
	DF <- object@normalizedTpmMatrix
	DF <- lapply(DF, Rle)
	DataFrame(DF)
})

#' @rdname CTSSnormalizedTpm

setMethod("CTSSnormalizedTpmDF",
signature(object = "CAGEexp"),
function (object){
  assays(object[["tagCountMatrix"]])$normalizedTpmMatrix
})

#' @name CTSSnormalizedTpmGR
#' @rdname CTSSnormalizedTpm
#' 
#' @param samples The name of sample(s) as reported by \code{sampleLabels(object)},
#'        or the number identifying the sample(s).
#'  
#' @export

setGeneric("CTSSnormalizedTpmGR", function(object, samples) {
  validSamples(object, samples)
  standardGeneric("CTSSnormalizedTpmGR")
  })

#' @rdname CTSSnormalizedTpm

setMethod( "CTSSnormalizedTpmGR", "CAGEr", function (object, samples) {
  gr <- CTSScoordinatesGR(object)
  score(gr) <- CTSSnormalizedTpmDF(object)[[samples]]
  gr <- gr[score(gr) != 0]
  .CTSS(gr)
})

#' @name CTSSclusteringMethod
#' 
#' @title Extracting CTSS clustering method from CAGEr objects.
#' 
#' @description Extracts the label of the method used for CTSS clustering into tag
#' clusters from \code{\link{CAGEr}} objects.
#' 
#' @param object A CAGEr object.
#' 
#' @return Returns a label of the method used for CTSS clustering.
#' 
#' @seealso \code{\link{clusterCTSS}}
#' 
#' @examples 
#' CTSSclusteringMethod(exampleCAGEset)
#' 
#' @author Vanja Haberle
#' @family CAGEr accessor methods
#' @export CTSSclusteringMethod

setGeneric("CTSSclusteringMethod", function(object) standardGeneric("CTSSclusteringMethod"))

#' @rdname CTSSclusteringMethod

setMethod("CTSSclusteringMethod", "CAGEset", function (object)
	object@clusteringMethod)

#' @rdname CTSSclusteringMethod

setMethod("CTSSclusteringMethod", "CAGEexp", function (object)
	metadata(object)$clusteringMethod)


#' @name tagClusters
#' 
#' @title Extract tag clusters (TCs) for individual CAGE experiments
#' 
#' @description Extracts tag clusters (TCs) produced by \code{\link{clusterCTSS}} function for
#' a specified CAGE experiment from a CAGEr object.
#' 
#' @param object A \code{\link{CAGEr}} object.
#' 
#' @param samples Label of the CAGE dataset (experiment, sample) for which to extract tag clusters.
#' If \code{samples = NULL}, a list of all the clusters for each sample is returned.
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
#' TC <- tagClusters(exampleCAGEset, "sample2", returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)
#' head(TC)
#' 
#' normalizeTagCount(exampleCAGEexp)
#' clusterCTSS( exampleCAGEexp, threshold = 50, thresholdIsTpm = TRUE
#'            , nrPassThreshold = 1, method = "distclu", maxDist = 20
#'            , removeSingletons = TRUE, keepSingletonsAbove = 100)
#' cumulativeCTSSdistribution(exampleCAGEexp, clusters = "tagClusters")
#' quantilePositions( exampleCAGEexp, clusters = "tagClusters"
#'                  , qLow = c(0.1,0.2), qUp = c(0.8,0.9))
#' TC <- tagClusters(exampleCAGEexp, "sample2", TRUE, 0.1, 0.9)
#' head(TC)

setGeneric( "tagClusters"
          , function( object, samples = NULL
                    , returnInterquantileWidth = FALSE, qLow = NULL, qUp = NULL) {
  validSamples(object, samples)
  standardGeneric("tagClusters")
})

.checkqLowUp <- function(qLow, qUp, object, sample) {
    if (is.null(qLow) | is.null(qUp))
      stop( "No quantiles specified! Please specify which quantile positions should be used to calculate width (qLow and qUp arguments)!")
    if (is.null(tagClustersQuantileLow(object)) & is.null(tagClustersQuantileUp(object)))
      stop("Interquantile width cannot be returned because no quantile positions have been calculated yet! Run 'quantilePositions()' first to get the positions of the desired quantiles!")
    if( (!(paste0("q_", qLow) %in% colnames(tagClustersQuantileLow(object)[[sample]]) &
           paste0("q_", qUp)  %in% colnames(tagClustersQuantileUp(object)[[sample]]))))
      stop("Interquantile width cannot be returned because specified quantile positions have not been calculated! Run 'quantilePositions()' again to get the positions of the desired quantiles!")
}
 
#' @rdname tagClusters

setMethod("tagClusters", "CAGEr", function (object, samples, returnInterquantileWidth, qLow, qUp){
  if (is.null(samples)) {
    tc.list <- lapply(sampleLabels(object), tagClusters, object = object, returnInterquantileWidth = returnInterquantileWidth, qLow = qLow, qUp = qUp)
    names(tc.list) <- sampleLabels(object)
    return(tc.list)
  }
  
  if (class(object) == "CAGEset") {
    tc <- object@tagClusters[[samples]]
  } else if (class(object) == "CAGEexp") {
    tc <- TCgranges2dataframe(metadata(object)$tagClusters[[samples]])
  } else {
    stop("Unsupported CAGEr class.")
  }
  
  if (returnInterquantileWidth) {
    .checkqLowUp(qLow, qUp, object, samples)
  	tc.w <- merge( tagClustersQuantileLow(object)[[samples]]
  	             , tagClustersQuantileUp(object)[[samples]])
  	tc.w <- tc.w[,c("cluster", paste0("q_", qLow), paste0("q_", qUp))]
  	tc.w$interquantile_width <- tc.w[,3] - tc.w[,2] + 1
  	tc <- merge(tc, tc.w)
  }
  tc
})

#' @rdname tagClusters
#' 
#' @param sample Label of one CAGE dataset (experiment, sample) for which to extract tag
#'        clusters. (For \code{tagClustersGR}, only one sample can be extracted.)
#' 
#' @export

setGeneric( "tagClustersGR"
          , function( object, sample = NULL
                    , returnInterquantileWidth = FALSE, qLow = NULL, qUp = NULL) {
  if (is.null(sample)) {
    tc.list <- GRangesList( lapply( sampleLabels(object)
                                  , tagClustersGR
                                  , object = object
                                  , returnInterquantileWidth = returnInterquantileWidth
                                  , qLow = qLow, qUp = qUp))
    names(tc.list) <- sampleLabels(object)
    return(tc.list)
  }
  validSamples(object, sample)
  standardGeneric("tagClustersGR")
})

#' @rdname tagClusters

setMethod( "tagClustersGR", "CAGEset"
         , function (object, sample, returnInterquantileWidth, qLow, qUp) {
  tc <- TCdataframe2granges(object@tagClusters[[sample]])
  if (returnInterquantileWidth) {
    .checkqLowUp(qLow, qUp, object, sample)
    qLowName <- paste0("q_", qLow)
    qUpName  <- paste0("q_", qUp)
    mcols(tc)[[qLowName]] <- tagClustersQuantileLow(object)[[sample]][,qLowName]
    mcols(tc)[[qUpName]]  <- tagClustersQuantileUp (object)[[sample]][,qUpName]
    tc$interquantile_width <- mcols(tc)[[qUpName]] - mcols(tc)[[qLowName]] + 1
  }
  .TagClusters(tc)
})

#' @rdname tagClusters

setMethod( "tagClustersGR", "CAGEexp"
         , function (object, sample, returnInterquantileWidth, qLow, qUp) {
  tc <- metadata(object)$tagClusters[[sample]]
  if (is.null(tc))
    stop( "No clusters found, run ", sQuote("clusterCTSS"), " first." , call. = FALSE)
  
  if (returnInterquantileWidth) {
    if (is.null(qLow) | is.null(qUp))
      stop( "No quantiles specified!  Set the ", sQuote("qLow")
          , " and ", sQuote("qUp"), "arguments.")
    qLowName <- paste0("q_", qLow)
    qUpName  <- paste0("q_", qUp)
    if(! all( c(qLowName, qUpName) %in% colnames(mcols(tc))))
      stop("Could not find quantile information.  Run ", sQuote("quantilePositions()"),  " first.")
    tc$interquantile_width <- mcols(tc)[[qUpName]] - mcols(tc)[[qLowName]] + 1
  }
  .TagClusters(tc)
})

#' @name filteredCTSSidx
#' @noRd

setGeneric("filteredCTSSidx", function(object) standardGeneric("filteredCTSSidx"))

setMethod("filteredCTSSidx", "CAGEset", function (object){
	object@filteredCTSSidx
})

setMethod("filteredCTSSidx", "CAGEexp", function (object){
  rowData(CTSStagCountSE(object))$filteredCTSSidx
})

#' @name tagClustersQuantile
#' @title Quantile metadata stored in CAGEr objects.
#' 
#' @description Accessor functions to quantile metadata.
#' 
#' @param object A \code{\link{CAGEr}} object.
#' @param samples Sample name(s), number(s) or \code{NULL} (default) for all samples.
#' @param value A list (one entry per sample) of data frames with multiple columns:
#'        \code{cluster} for the cluster ID, and then \code{q_0.n} where \code{0.n}
#'        indicates a quantile.
NULL

#' @name tagClustersQuantileLow
#' @rdname tagClustersQuantile

setGeneric("tagClustersQuantileLow", function(object, samples = NULL) {
  validSamples(object, samples)
  standardGeneric("tagClustersQuantileLow")
})

#' @rdname tagClustersQuantile

setMethod("tagClustersQuantileLow", "CAGEset", function (object, samples) {
	if (is.null(samples)) return(object@tagClustersQuantileLow)
  object@tagClustersQuantileLow[[samples]]
})

#' @rdname tagClustersQuantile

setMethod("tagClustersQuantileLow", "CAGEexp", function (object, samples) {
  stop("Not supported for CAGEexp; get quantile data via ", sQuote("tagClustersGR") , ".")
})


#' @name tagClustersQuantileUp
#' @rdname tagClustersQuantile

setGeneric("tagClustersQuantileUp", function(object, samples = NULL) {
  validSamples(object, samples)
  standardGeneric("tagClustersQuantileUp")
})

#' @rdname tagClustersQuantile

setMethod("tagClustersQuantileUp", "CAGEset", function (object, samples) {
  if (is.null(samples)) return(object@tagClustersQuantileUp)
  object@tagClustersQuantileUp[[samples]]
})

#' @rdname tagClustersQuantile

setMethod("tagClustersQuantileUp", "CAGEexp", function (object, samples) {
  stop("Not supported for CAGEexp; get quantile data via ", sQuote("tagClustersGR") , ".")
})


#' @name getConsensusClusters
#' @noRd

setGeneric("getConsensusClusters", function(object) standardGeneric("getConsensusClusters"))

setMethod("getConsensusClusters", "CAGEset", function (object){
	object@consensusClusters
})

setMethod("getConsensusClusters", "CAGEexp", function (object){
  cc <- rowRanges(consensusClustersSE(object))
  cc$consensus.cluster <- names(cc)
  CCgranges2dataframe(cc)
})


#' @name consensusClustersGR
#' @rdname consensusClusters
#' @export

setGeneric("consensusClustersGR", function(object, samples = NULL) {
  validSamples(object, samples)
  standardGeneric("consensusClustersGR")})

#' @rdname consensusClusters

setMethod("consensusClustersGR", "CAGEset", function (object, samples) {
  if (is.null(samples)) return(CCdataframe2granges(object@consensusClusters))
  CCdataframe2granges(object@consensusClusters)})

#' @rdname consensusClusters

setMethod("consensusClustersGR", "CAGEexp", function (object, samples) {
  if(is.null(experiments(object)$consensusClusters))
    stop("No consensus clusters found.  See ", sQuote("?aggregateTagClusters"), " on how to create them.")
  .ConsensusClusters(rowRanges(consensusClustersSE(object)))
})


#' @name consensusClustersSE
#' @rdname consensusClusters
#' @export

setGeneric("consensusClustersSE", function(object) standardGeneric("consensusClustersSE"))

#' @rdname consensusClusters

setMethod("consensusClustersSE", "CAGEset", function (object)
	stop("CAGEset objects not supported"))

#' @rdname consensusClusters

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
#' exampleCAGEexp$group <- c("a", "a", "b", "b", "a")
#' consensusClustersDESeq2(exampleCAGEexp, ~group)
#' 
#' @export

setGeneric( "consensusClustersDESeq2"
          , function(object, design) standardGeneric("consensusClustersDESeq2"))

#' @rdname consensusClustersDESeq2

setMethod( "consensusClustersDESeq2", "CAGEset"
         , function (object, design) stop("Not implemented for the CAGEset class."))

#' @rdname consensusClustersDESeq2

setMethod( "consensusClustersDESeq2", "CAGEexp"
         , function (object, design) {
  if (! requireNamespace("DESeq2"))
    stop("This function requires the ", dQuote("DESeq2"), " package; please install it.")
  DESeq2::DESeqDataSetFromMatrix( countData = assay(consensusClustersSE(object))
                                , colData   = colData(object)
                                , rowData   = consensusClustersGR(object)
                                , design    = design)
})

#' @name consensusClusters
#' 
#' @title Get or set consensus clusters from CAGEr objects
#' 
#' @description Extracts the information on consensus clusters from a CAGEr object.
#' 
#' @param object A \code{\link{CAGEr}} object.
#' 
#' @param samples Optional. Label of the CAGE dataset (experiment, samples) for which to extract
#' sample-specific information on consensus clusters.
#' 
#' @param returnInterquantileWidth Should the interquantile width of consensus clusters in
#' specified sample be returned.  Used only when \code{samples} argument is specified, otherwise
#' ignored.
#' 
#' @param qLow Position of which quantile should be used as a left (lower) boundary when
#' calculating interquantile width.  Used only when \code{samples} argument is specified and
#' \code{returnInterquantileWidth = TRUE}, otherwise ignored.
#' 
#' @param qUp Position of which quantile should be used as a right (upper) boundary when
#' calculating interquantile width.  Used only when \code{samples} argument is specified and
#' \code{returnInterquantileWidth = TRUE}, otherwise ignored.
#' 
#' @return Returns a \code{data.frame} with information on consensus clusters, including genomic
#' coordinates.  When \code{samples} argument is NOT specified, total CAGE signal across all CAGE
#' datasets (samples) is returned in the \code{tpm} column.  When \code{samples} argument is
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
#' head(consensusClusters(exampleCAGEset))
#' head(consensusClusters(exampleCAGEset, samples = "sample2"))
#' 
#' clusterCTSS(exampleCAGEexp)
#' aggregateTagClusters(exampleCAGEexp)
#' head(consensusClusters(exampleCAGEexp))
#' consensusClustersGR(exampleCAGEexp, 2)
#' 
#' @export

setGeneric(
name="consensusClusters",
def=function(object, samples=NULL, returnInterquantileWidth = FALSE, qLow = NULL, qUp = NULL){
    standardGeneric("consensusClusters")
})

#' @rdname consensusClusters

setMethod("consensusClusters",
signature(object = "CAGEr"),
function (object, samples = NULL, returnInterquantileWidth = FALSE, qLow = NULL, qUp = NULL){
	
    if(is.null(samples)){
        return(getConsensusClusters(object))
    }else if(samples %in% sampleLabels(object)){
      if(class(object) == "CAGEexp")
        stop( sQuote("samples"), " option not supported for CAGEexp objects.  Use "
            , sQuote("consensusClustersGR()"), " instead.")
        
        cc.s <- cbind(cluster = as.integer(rownames(consensusClustersTpm(object))), tpm = consensusClustersTpm(object)[,samples])
        
        if(returnInterquantileWidth & (length(qLow) == 0 | length(qUp) == 0)){
            stop("No quantiles specified! Please specify which quantile positions should be used to calculate width (qLow and qUp arguments)!")
        }else if(returnInterquantileWidth & (length(consensusClustersQuantileLow(object))==0 & length(consensusClustersQuantileUp(object))==0)){
            stop("Interquantile width cannot be returned because no quantile positions for consensus clusters have been calculated yet! Run 'quantilePositions()' first to get the positions of the desired quantiles!")
        }else if(returnInterquantileWidth & (!(paste("q_", qLow, sep = "") %in% colnames(consensusClustersQuantileLow(object)[[samples]]) & paste("q_", qUp, sep = "") %in% colnames(consensusClustersQuantileUp(object)[[samples]])))){
            stop("Interquantile width cannot be returned because specified quantile positions have not been calculated for consensus clusters! Run 'quantilePositions()' again to get the positions of the desired quantiles!")
        }else if(returnInterquantileWidth){
            cc <- getConsensusClusters(object)

            cc.cumsum <- CTSScumulativesCC(object)[[samples]]
            a <- lapply(cc.cumsum, function(x) {.get.dominant.ctss(as.numeric(x), isCumulative = T)})
            b <- data.frame(consensus.cluster = as.integer(names(a)), dominant_ctss = unlist(a))
            #cc <- merge(b, cc.s, by.x = 1, by.y = 1, all.x = T, all.y = F)
            
            cc <- merge(cc[,-which(colnames(cc) == "tpm")], b, by.x = 1, by.y = 1, all.x = F, all.y = T)
            cc$dominant_ctss <- cc$start + cc$dominant_ctss
            
            cc <- merge(cc, cc.s, by.x = 1, by.y = 1, all.x = T, all.y = F)
            
            ctss <- CTSSnormalizedTpm(object)[,c("chr", "pos", "strand", samples)]
            cc <- merge(cc, ctss, by.x = c("chr", "strand", "dominant_ctss"), by.y = c("chr", "strand", "pos"), all.x = T, all.y = F)
            colnames(cc)[ncol(cc)] <- "tpm.dominant_ctss"
            cc <- cc[,c("consensus.cluster", "chr", "start", "end", "strand", "dominant_ctss", "tpm", "tpm.dominant_ctss")]
            
            cc.w <- merge(consensusClustersQuantileLow(object)[[samples]], consensusClustersQuantileUp(object)[[samples]])
            cc.w <- cc.w[,c(1, which(colnames(cc.w) == paste("q_", qLow, sep = "")), which(colnames(cc.w) == paste("q_", qUp, sep = "")))]
            cc.w$interquantile_width <- cc.w[,3] - cc.w[,2] + 1
            cc <- merge(cc, cc.w, by.x = "consensus.cluster", by.y = "cluster", all.x = T)
        }else{
            cc <- getConsensusClusters(object)
            cc <- merge(cc[,-which(colnames(cc) == "tpm")], cc.s, by.x = 1, by.y = 1, all.x = T, all.y = F)
            cc <- subset(cc, cc$tpm>0)
        }
        return(cc)
    }else{
        stop("Provided 'samples' not in the CAGE set! Check sampleLabels()")
    }

})

#' @name consensusClustersQuantile
#' @title Quantile metadata stored in CAGEr objects.
#' 
#' @description Accessors for consensus cluster quantile data in CAGEr objects.
#' 
#' @param object A \code{\link{CAGEr}} object.
#' @param value A list (one entry per sample) of data frames with multiple columns:
#'        \code{cluster} for the cluster ID, and then \code{q_0.n} where \code{0.n}
#'        indicates a quantile.
#' @param samples Sample name(s), number(s) or \code{NULL} (default) for all samples.
NULL

#' @name consensusClustersQuantileLow
#' @rdname consensusClustersQuantile

setGeneric("consensusClustersQuantileLow", function(object, samples = NULL) {
  validSamples(object, samples)
  standardGeneric("consensusClustersQuantileLow")
})

#' @rdname consensusClustersQuantile

setMethod("consensusClustersQuantileLow", "CAGEset", function (object, samples) {
	if (is.null(samples)) return(object@consensusClustersQuantileLow)
  object@consensusClustersQuantileLow[[samples]]
})

#' @rdname consensusClustersQuantile

setMethod("consensusClustersQuantileLow", "CAGEexp", function (object, samples) {
  if (is.null(samples)) return(metadata(object)$consensusClustersQuantileLow)
  metadata(object)$consensusClustersQuantileLow[[samples]]
})


#' @name consensusClustersQuantileUp
#' @rdname consensusClustersQuantile

setGeneric("consensusClustersQuantileUp", function(object, samples = NULL) {
  validSamples(object, samples)
  standardGeneric("consensusClustersQuantileUp")
})

#' @rdname consensusClustersQuantile

setMethod("consensusClustersQuantileUp", "CAGEset", function (object, samples) {
  if (is.null(samples)) return(object@consensusClustersQuantileUp)
  object@consensusClustersQuantileUp[[samples]]
})

#' @rdname consensusClustersQuantile

setMethod("consensusClustersQuantileUp", "CAGEexp", function (object, samples) {
  if (is.null(samples)) return(metadata(object)$consensusClustersQuantileUp)
  metadata(object)$consensusClustersQuantileUp[[samples]]
})



#' @name CTSScumulativesTagClusters
#'  
#' @title Get/set CTSS cumulative TC or CC data
#' 
#' @description Accessor function.
#' 
#' @param object A \code{\link{CAGEset}} or \code{\link{CAGEset}} object.
#' @param samples One or more valid sample names.
#' 
#' @return List of numeric Rle.
#' 
#' @family CAGEr clusters functions
#' @family CAGEr accessor methods
#' 
#' @export

setGeneric("CTSScumulativesTagClusters", function(object, samples = NULL) {
  validSamples(object, samples)
  standardGeneric("CTSScumulativesTagClusters")
})

#' @rdname CTSScumulativesTagClusters

setMethod("CTSScumulativesTagClusters", "CAGEset", function (object, samples){
    if (is.null(samples)) return(object@CTSScumulativesTagClusters)
	object@CTSScumulativesTagClusters[[samples]]
})

#' @rdname CTSScumulativesTagClusters

setMethod("CTSScumulativesTagClusters", "CAGEexp", function (object, samples){
  tc <- metadata(object)$CTSScumulativesTagClusters
  if (is.null(tc))
    stop( "No cumulative sums found, run ", sQuote("cumulativeCTSSdistribution"), " first."
        , call. = FALSE)
  if (is.null(samples)) return(tc)
  tc[[samples]]
})

#' @rdname CTSScumulativesTagClusters

setGeneric("CTSScumulativesCC", function(object, samples = NULL) {
  validSamples(object, samples)
  standardGeneric("CTSScumulativesCC")
})

#' @rdname CTSScumulativesTagClusters

setMethod("CTSScumulativesCC", "CAGEset", function (object, samples) {
  if (is.null(samples)) return(object@CTSScumulativesConsensusClusters)
	object@CTSScumulativesConsensusClusters[[samples]]
})

#' @rdname CTSScumulativesTagClusters

setMethod("CTSScumulativesCC", "CAGEexp", function (object, samples) {
  cc <- metadata(object)$CTSScumulativesConsensusClusters
  if (is.null(cc))
    stop( "No cumulative sums found, run ", sQuote("cumulativeCTSSdistribution"), " first."
        , call. = FALSE)
  if (is.null(samples)) return(cc)
  cc[[samples]]
})

#' @name consensusClustersTpm
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
#' clusters.tpm <- consensusClustersTpm(exampleCAGEset)
#' head(clusters.tpm)
#' consensusClustersTpmDF(exampleCAGEset)
#' 
#' @importFrom SummarizedExperiment assays
#' 
#' @export consensusClustersTpm

setGeneric("consensusClustersTpm", function(object) standardGeneric("consensusClustersTpm"))

#' @rdname consensusClustersTpm

setMethod("consensusClustersTpm",
signature(object = "CAGEset"),
function (object){
	object@consensusClustersTpmMatrix
})

#' @rdname consensusClustersTpm

setMethod("consensusClustersTpm",
signature(object = "CAGEexp"),
function (object){
  consensusClustersSE <- object[["consensusClusters"]]
  as.matrix(assays(consensusClustersSE)[["normalized"]])
})


#' @name consensusClustersTpmDf
#' @rdname consensusClustersTpm
#' @export

setGeneric("consensusClustersTpmDf", function(object) standardGeneric("consensusClustersTpmDf"))

#' @rdname consensusClustersTpm

setMethod("consensusClustersTpmDf", "CAGEset", function (object)
	as.data.frame(object@consensusClustersTpmMatrix))

#' @rdname consensusClustersTpm

setMethod("consensusClustersTpmDf", "CAGEexp", function (object) {
  consensusClustersSE <- object[["consensusClusters"]]
  as.data.frame(assays(consensusClustersSE)[["normalized"]])
})


#' @name consensusClustersTpmDF
#' @rdname consensusClustersTpm
#' @export

setGeneric("consensusClustersTpmDF", function(object) standardGeneric("consensusClustersTpmDF"))

#' @rdname consensusClustersTpm

setMethod("consensusClustersTpmDF", "CAGEset", function (object)
	DataFrame(object@consensusClustersTpmMatrix))

#' @rdname consensusClustersTpm

setMethod("consensusClustersTpmDF", "CAGEexp", function (object) {
  consensusClustersSE <- object[["consensusClusters"]]
  assays(consensusClustersSE)[["normalized"]]
})


#' @name expressionClasses
#' 
#' @title Extract labels of expression classes
#' 
#' @description Retrieves labels of expression classes of either individual CTSSs or consensus
#' clusters from a CAGEset object.
#' 
#' @param object A \code{\link{CAGEset}} object.
#' 
#' @param what Which level of expression clustering should be used. Can be either
#'        \code{"CTSS"} to extract labels of expression classes of individual CTSSs or
#'        \code{"consensusClusters"} to extract labels of expression classes of consensus
#'        clusters.
#' 
#' @return Returns character vector of labels of expression classes.  The number of labels
#' matches the number of expression clusters returned by \code{\link{getExpressionProfiles}}
#' function.
#' 
#' @seealso \code{\link{getExpressionProfiles}} \code{\link{plotExpressionProfiles}}
#'          \code{\link{extractExpressionClass}}
#'          
#' @examples
#' 
#' exprClasses <- expressionClasses(exampleCAGEset, what = "CTSS")
#' exprClasses
#' 
#' @export

setGeneric( "expressionClasses"
          , function(object, what = c("CTSS", "consensusClusters"))
              standardGeneric("expressionClasses"))

#' @rdname expressionClasses

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

#' @rdname expressionClasses

setMethod("expressionClasses", "CAGEexp", function (object, what) {
  switch (match.arg(what),
    CTSS = {
      classes <- CTSScoordinatesGR(object)$expressionClasses
      if(is.null(classes)) stop("No expression clustering of CTSSs has been done yet!")},
    consensusClusters = {
      classes <- consensusClustersGR(object)$expressionClasses
      if(is.null(classes)) stop("No expression clustering of consensus clusters has been done yet!")})
  unique(sort(classes))
})


#' @name CTSSexpressionClusteringMethod
#' @noRd 

setGeneric("CTSSexpressionClusteringMethod", function(object) 
  standardGeneric("CTSSexpressionClusteringMethod"))

setMethod("CTSSexpressionClusteringMethod", "CAGEset", function (object)
  object@CTSSexpressionClusteringMethod)

setMethod("CTSSexpressionClusteringMethod", "CAGEexp", function (object)
  metadata(object)$CTSSexpressionClusteringMethod)



#' @name CTSSexpressionClasses
#' @noRd

setGeneric("CTSSexpressionClasses", function(object) standardGeneric("CTSSexpressionClasses"))

setMethod("CTSSexpressionClasses", "CAGEset", function (object)
  object@CTSSexpressionClasses)

setMethod("CTSSexpressionClasses", "CAGEexp", function (object)
  metadata(object)$CTSSexpressionClasses)


#' @name consensusClustersExpressionClasses
#' @noRd 

setGeneric("consensusClustersExpressionClasses", function(object)
  standardGeneric("consensusClustersExpressionClasses"))

setMethod("consensusClustersExpressionClasses", "CAGEset", function (object)
  object@consensusClustersExpressionClasses)

setMethod("consensusClustersExpressionClasses", "CAGEexp", function (object)
  metadata(object)$consensusClustersExpressionClasses)


#' @name consensusClustersExpressionClusteringMethod
#' @noRd 

setGeneric("consensusClustersExpressionClusteringMethod", function(object) 
  standardGeneric("consensusClustersExpressionClusteringMethod"))

setMethod("consensusClustersExpressionClusteringMethod", "CAGEset", function (object)
  object@consensusClustersExpressionClusteringMethod)

setMethod("consensusClustersExpressionClusteringMethod", "CAGEexp", function (object)
  metadata(object)$consensusClustersExpressionClusteringMethod)


#' @name GeneExpSE
#' 
#' @title Retreives the SummarizedExperiment containing gene expression levels.
#' 
#' @description Get or set a \code{SummarizedExperiment} using the gene expression
#' data in the experiment slot \code{geneExpMatrix} and the sample metadata
#' of the \code{\link{CAGEexp}} object.
#' 
#' @param object A \code{\link{CAGEexp}} object.
#' 
#' @family CAGEr accessor methods
#' 
#' @author Charles Plessy
#' 
#' @examples 
#' GeneExpSE(exampleCAGEexp)
#' 
#' @export

setGeneric("GeneExpSE", function(object) standardGeneric("GeneExpSE"))

#' @rdname GeneExpSE

setMethod("GeneExpSE", "CAGEset", function (object)
	stop("Not implemented for the CAGEset class."))

#' @rdname GeneExpSE

setMethod("GeneExpSE", "CAGEexp", function (object)
  experiments(object)$geneExpMatrix)


#' @name GeneExpDESeq2
#' 
#' @title Export gene expression data for DESeq2 analysis
#' 
#' @description Creates a \code{DESeqDataSet} using the gene expression
#' data in the experiment slot \code{geneExpMatrix} and the sample metadata
#' of the \code{\link{CAGEexp}} object.  The formula must be built using factors
#' already present in the sample metadata.
#' 
#' @param object A \code{\link{CAGEexp}} object.
#' @param design A formula for the DESeq2 analysis.
#' 
#' @author Charles Plessy
#' 
#' @seealso \code{DESeqDataSet} in the \code{DESeq2} package.
#' 
#' @family CAGEr gene expression analysis functions
#' @family CAGEr accessor methods
#' 
#' @examples
#' exampleCAGEexp$group <- factor(c("a", "a", "b", "b", "a"))
#' GeneExpDESeq2(exampleCAGEexp, ~group)
#' 
#' @export

setGeneric( "GeneExpDESeq2"
          , function(object, design) standardGeneric("GeneExpDESeq2"))

#' @rdname GeneExpDESeq2

setMethod( "GeneExpDESeq2", "CAGEset"
         , function (object, design) stop("Not implemented for the CAGEset class."))

#' @rdname GeneExpDESeq2

setMethod( "GeneExpDESeq2", "CAGEexp"
         , function (object, design) {
  if (! requireNamespace("DESeq2"))
    stop("This function requires the ", dQuote("DESeq2"), " package; please install it.")
  DESeq2::DESeqDataSetFromMatrix( countData = assay(GeneExpSE(object))
                                , colData   = colData(object)
                                , rowData   = rowData(GeneExpSE(object))
                                , design    = design)
})