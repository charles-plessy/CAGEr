#' @include AggregationFunctions.R ClusteringFunctions.R CAGEr.R CTSS.R

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
#' @details \code{\link{CAGEexp}} objects constructed with `NULL` in place
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

setGeneric("genomeName", function(object) standardGeneric("genomeName"))

#' @rdname genomeName

setMethod("genomeName", "CAGEset", function (object)	object@genomeName)

#' @rdname genomeName

setMethod("genomeName", "CAGEexp", function (object) metadata(object)$genomeName)

#' @rdname genomeName

setMethod("genomeName", "CTSS", function (object) metadata(object)$genomeName)


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
#' @title Genomic coordinates of TSSs from a `CAGEr` object
#' 
#' @description Extracts the genomic coordinates of all detected TSSs
#' from [CAGEset] and [CAGEexp] objects.
#' 
#' @param object A `CAGEset` or `CAGEexp` object.
#' 
#' @return [CTSScoordinates()] returns a `data.frame` with genomic coordinates
#' of all TSSs.  The `pos` column contains 1-based coordinate of the TSS.
#' 
#' @seealso
#' [`getCTSS`]
#' 
#' @examples
#' CTSS <- CTSScoordinates(exampleCAGEset)
#' head(CTSS)
#'
#' @author Vanja Haberle
#' @family CAGEr accessor methods
#' @export

setGeneric("CTSScoordinates", function(object) standardGeneric("CTSScoordinates"))

#' @rdname CTSScoordinates

setMethod("CTSScoordinates", "CAGEset", function (object)
  object@CTSScoordinates)

#' @rdname CTSScoordinates

setMethod("CTSScoordinates", "CAGEexp", function (object) {
  gr <- rowRanges(experiments(object)$tagCountMatrix)
  data.frame( chr = as.character(seqnames(gr))
            , pos = start(gr)
            , strand = as.character(strand(gr))
            , stringsAsFactors = FALSE)
})


#' @name CTSScoordinatesGR
#' @rdname CTSScoordinates
#' 
#' @return `CTSScoordinatesGR` returns the coordinates as a [CTSS()] object
#' wrapping genomic ranges.  A `filteredCTSSidx` column metadata will be present
#' if [clusterCTSS()] was ran earlier.
#' 
#' @seealso [`clusterCTSS`]
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

setGeneric("CTSScoordinatesGR", function(object) standardGeneric("CTSScoordinatesGR"))

#' @rdname CTSScoordinates

setMethod( "CTSScoordinatesGR", "CAGEset", function (object){
  ctssCoord <- object@CTSScoordinates
  ctssCoord <- GPos(stitch = FALSE
                   ,GRanges( ctssCoord$chr
                           , IRanges(ctssCoord$pos, width = 1)
                           , ctssCoord$strand))
  if(!identical(object@filteredCTSSidx, logical()))
    ctssCoord$filteredCTSSidx <- object@filteredCTSSidx
  if(!identical(CTSSexpressionClasses(object), character()))
    if (length(CTSSexpressionClasses(object)) == length(ctssCoord)) {
      ctssCoord$CTSSexpressionClasses <- CTSSexpressionClasses(object)
    } else {
      warning("Skipping expression classes: not same length as CTSS.")
    }
  .CTSS(ctssCoord, bsgenomeName = genomeName(object))
})

#' @rdname CTSScoordinates

setMethod("CTSScoordinatesGR", "CAGEexp", function (object)
  rowRanges(CTSStagCountSE(object)))


#' @name CTSStagCount
#' 
#' @title Extract raw CAGE TSSs expression tables from [`CAGEr`] objects
#' 
#' @description Extracts the tag count for all detected TSSs in all CAGE datasets
#'              from [`CAGEset`] and [`CAGEexp`] objects.
#' 
#' @param object A `CAGEr` object.
#' @param samples (For `CTSStagCountGR`.) Name(s) or number(s) identifying sample(s).
#'  
#' @return Returns an object with number of CAGE tags supporting each TSS
#' (rows) in every CAGE dataset (columns).  The class of the object depends on the
#' function being called:
#' 
#' * `CTSStagCount`:   A [`data.frame`], containing coordinates and expression values.
#' * `CTSStagCountDf`: A `data.frame``, containing only expression values.
#' * `CTSStagCountDF`: A [`DataFrame`] of [`Rle`] integers.
#' * `CTSStagCountDA`: A [`DelayedArray`] wrapping a `DataFrame` of `Rle` integers.
#' * `CTSStagCountSE`: A [`RangedSummarizedExperiment`]` containing a `DataFrame`
#'    of `Rle` integers.
#' * `CTSStagCountGR`: A `CTSS` object (wrapping `GRanges`) containing a `score`
#'     column indicating expression values for a given sample.
#' 
#' @seealso [getCTSS()]
#' 
#' @examples
#' head(CTSStagCount(exampleCAGEset))
#' head(CTSStagCountDf(exampleCAGEset))
#' CTSStagCountDF(exampleCAGEset)
#' CTSStagCountDA(exampleCAGEset)
#' 
#' head(CTSStagCount(exampleCAGEexp))
#' head(CTSStagCountDf(exampleCAGEexp))
#' CTSStagCountDF(exampleCAGEset)
#' CTSStagCountDA(exampleCAGEset)
#' 
#' @author Vanja Haberle
#' @author Charles Plessy
#' 
#' @family CAGEr accessor methods
#' @family CAGEr CTSS methods
#' @export

setGeneric("CTSStagCount", function(object) standardGeneric("CTSStagCount"))

#' @rdname CTSStagCount

setMethod("CTSStagCount", "CAGEset", function (object)
	cbind(object@CTSScoordinates, object@tagCountMatrix))

#' @rdname CTSStagCount

setMethod( "CTSStagCount", "CAGEexp", function (object)
  cbind( CTSScoordinates(object)
       , as.data.frame(lapply(assay(experiments(object)$tagCountMatrix), as.integer))))


#' @name CTSStagCountDf
#' @rdname CTSStagCount
#' 
#' @export

setGeneric("CTSStagCountDf", function(object) standardGeneric("CTSStagCountDf"))

#' @rdname CTSStagCount

setMethod("CTSStagCountDf", "CAGEset", function (object) 	object@tagCountMatrix)

#' @rdname CTSStagCount

setMethod("CTSStagCountDf", "CAGEexp", function (object) 
  as.data.frame(lapply(assay(experiments(object)$tagCountMatrix), as.integer)))


#' @name CTSStagCountDF
#' @rdname CTSStagCount
#'  
#' @export

setGeneric("CTSStagCountDF", function(object) standardGeneric("CTSStagCountDF"))

#' @rdname CTSStagCount

setMethod("CTSStagCountDF", "CAGEset", function (object) {
	DF <- object@tagCountMatrix
	DF <- lapply(DF, as.integer)
	DF <- lapply(DF, Rle)
	DataFrame(DF)
})

#' @rdname CTSStagCount

setMethod("CTSStagCountDF", "CAGEexp", function (object)
  assay(CTSStagCountSE(object)))


#' @name CTSStagCountDA
#' @rdname CTSStagCount
#' 
#' @import DelayedArray DelayedArray
#' @export

setGeneric("CTSStagCountDA", function(object) tandardGeneric("CTSStagCountDA"))

#' @rdname CTSStagCount

setMethod("CTSStagCountDA", signature(object = "CAGEr"), function (object)
  DelayedArray(CTSStagCountDF(object)))


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
  if (is.character(samples)) samples <- which(sampleLabels(object) == samples)
  gr <- CTSScoordinatesGR(object)
  score(gr) <- CTSStagCountDF(object)[[samples]]
  gr <- gr[score(gr) != 0]
  gr <- .CTSS(gr, bsgenomeName = genomeName(object))
  sampleLabels(gr) <- sampleLabels(object)[samples]
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
  .CTSS(gr, bsgenomeName = genomeName(object))
})

#' @name CTSSclusteringMethod
#' 
#' @title Get /set CTSS clustering method
#' 
#' @description Returns or sets the name of the method that was used make tag
#' clusters from the CTSSs of a \code{\link{CAGEr}} object.
#' 
#' @param object A CAGEr object.
#' 
#' @seealso \code{\link{clusterCTSS}}
#' @family CAGEr accessor methods
#' @family CAGEr clusters functions
#' 
#' @author Vanja Haberle
#' @author Charles Plessy
#' 
#' @examples 
#' CTSSclusteringMethod(exampleCAGEset)
#' CTSSclusteringMethod(exampleCAGEexp)
#' 
#' @export CTSSclusteringMethod

setGeneric("CTSSclusteringMethod", function(object) standardGeneric("CTSSclusteringMethod"))

#' @rdname CTSSclusteringMethod

setMethod("CTSSclusteringMethod", "CAGEset", function (object)
	object@clusteringMethod)

#' @rdname CTSSclusteringMethod

setMethod("CTSSclusteringMethod", "GRangesList", function (object)
	metadata(object)$clusteringMethod)

#' @rdname CTSSclusteringMethod

setMethod("CTSSclusteringMethod", "CAGEexp", function (object)
  CTSSclusteringMethod(metadata(object)$tagClusters))
  # extrat directly TCs from metadata slot because tagClustersGR does more that
  # is not needed here.


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
#' head(tagClusters( exampleCAGEset, "sample2"
#'                 , returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9))

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
    tc.list <- lapply(sampleList(object), tagClusters, object = object, returnInterquantileWidth = returnInterquantileWidth, qLow = qLow, qUp = qUp)
    return(tc.list)
  }
  
  if (inherits(object, "CAGEset")) {
    tc <- object@tagClusters[[samples]]
  } else if (inherits(object, "CAGEexp")) {
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
#' @examples
#' tagClustersGR(exampleCAGEexp, "Zf.high", TRUE, 0.1, 0.9)
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
    metadata(tc.list)$clusteringMethod <- CTSSclusteringMethod(object)
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
    mcols(tc)[[qLowName]] <- tagClustersQuantileLow(object)[[sample]][,qLowName] - start(tc)
    mcols(tc)[[qUpName]]  <- tagClustersQuantileUp (object)[[sample]][,qUpName]  - start(tc)
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
#' @param q A single quantile (not a list)
#' @param value A list (one entry per sample) of data frames with multiple columns:
#'        \code{cluster} for the cluster ID, and then \code{q_0.n} where \code{0.n}
#'        indicates a quantile.
#' 
#' @return Returns a \code{data.frame} where the first column gives cluster
#' names and the next columns give quantile positions, in \emph{zero-based}
#' chromosome coordinates (because the tag clusters in CAGEset objects are
#' represented in zero-based coordinates as well)).

setGeneric("tagClustersQuantile", function(object, samples = NULL, q = NULL)
  standardGeneric("tagClustersQuantile"))

.getClustersQuantile <- function (object, q) {
  qName <- paste0("q_", q)
  if (! all(qName %in% colnames(mcols(object))))
    stop( "At least one of the quantiles "
        , paste(sQuote(qName), collapse=" or ")
        , " was not found.")
  tcq <- mcols(object)[, qName, drop = FALSE]
	tcq <- data.frame(lapply(tcq, decode))
	tcq <- tcq + start(object) - 1
	cbind(cluster = names(object), tcq)
}

#' @rdname tagClustersQuantile

setMethod("tagClustersQuantile", "TagClusters", function (object, samples, q) {
  if (is.null(q))
    stop("Indicate which quantile(s) to extract from this TagClusters object.")
  if (! is.null(samples))
    stop(sQuote("samples"), " must be NULL.")
  .getClustersQuantile(object, q)
})

#' @rdname tagClustersQuantile

setMethod("tagClustersQuantile", "CAGEexp", function (object, samples, q) {
  validSamples(object, samples)
  if (length(samples) > 1)
    stop("Multiple samples not supported for CAGEexp objects(all or just one).")
  if (is.null(q))
    stop("A single quantile must be chosen for CAGEexp objects.")
  if (is.null(samples)) {
    lapply( sampleList(object)
          , tagClustersQuantile
          , object = object
          , q = q)
  } else {
    tagClustersQuantile(tagClustersGR(object, samples), q = q)
  }
})

#' @name tagClustersQuantileLow
#' @rdname tagClustersQuantile

setGeneric("tagClustersQuantileLow", function(object, samples = NULL, q = NULL) {
  validSamples(object, samples)
  standardGeneric("tagClustersQuantileLow")
})

#' @rdname tagClustersQuantile

setMethod("tagClustersQuantileLow", "CAGEset", function (object, samples, q) {
  if (is.null(object@tagClustersQuantileLow))
    stop( "No data for given quantile positions! "
        , "Run 'quantilePositions()' function for desired quantiles first!")
  res <- object@tagClustersQuantileLow
  if (is.null(samples)) {
    if (!is.null(q) & !all(sapply(res, function(x, q) paste0("q_", q) %in% names(x), q = q)))
      stop( "Low quantile not found! "
          , "Run 'quantilePositions()' function for desired quantiles first!")
    res
  } else {
    res <- res[[samples]]
    if (!is.null(q) & !paste0("q_", q) %in% names(res))
      stop( "Low quantile not found! "
          , "Run 'quantilePositions()' function for desired quantiles first!")
    res
  }
})

#' @rdname tagClustersQuantile

setMethod("tagClustersQuantileLow", "CAGEexp", function (object, samples, q)
  tagClustersQuantile(object = object, samples = samples, q = q))


#' @name tagClustersQuantileUp
#' @rdname tagClustersQuantile

setGeneric("tagClustersQuantileUp", function(object, samples = NULL, q = NULL) {
  validSamples(object, samples)
  standardGeneric("tagClustersQuantileUp")
})

#' @rdname tagClustersQuantile

setMethod("tagClustersQuantileUp", "CAGEset", function (object, samples, q) {
  if (is.null(object@tagClustersQuantileUp))
    stop( "No data for given quantile positions! "
        , "Run 'quantilePositions()' function for desired quantiles first!")
  res <- object@tagClustersQuantileUp
  if (is.null(samples)) {
    if (!is.null(q) & !all(sapply(res, function(x, q) paste0("q_", q) %in% names(x), q = q)))
      stop( "Low quantile not found! "
          , "Run 'quantilePositions()' function for desired quantiles first!")
    res
  } else {
    res <- res[[samples]]
    if (!is.null(q) & !paste0("q_", q) %in% names(res))
      stop( "Low quantile not found! "
          , "Run 'quantilePositions()' function for desired quantiles first!")
    res
  }
})

#' @rdname tagClustersQuantile

setMethod("tagClustersQuantileUp", "CAGEexp", function (object, samples, q)
  tagClustersQuantile(object = object, samples = samples, q = q))


#' @name consensusClustersGR
#' @rdname consensusClusters
#' @return `consensusClustersGR` returns a [`ConsensusClusters`] object, which
#' is a [`GRanges`] wrapper containing similar information as the data frame
#' returned by `consensusClustersGR`.  The `score` columns indicates the
#' normalised expression value of each cluster, either across all samples
#' (`sample = NULL`), or for the selected sample.  The `tpm` column provides
#' the same information for compatibility with `CAGEset` objects but may be
#' removed in the future.
#' 
#' @export

setGeneric( "consensusClustersGR"
          , function( object
                    , sample = NULL
                    , returnInterquantileWidth = FALSE
                    , qLow = NULL, qUp = NULL) {
  validSamples(object, sample)
  standardGeneric("consensusClustersGR")})

#' @rdname consensusClusters

setMethod( "consensusClustersGR", "CAGEset"
         , function (object, sample, returnInterquantileWidth, qLow, qUp) {
           
  if (returnInterquantileWidth)
    stop( dQuote("returnInterquantileWidth")
        ," argument is not supported for CAGEset objects.")
  if (!is.null(qLow))
    stop(dQuote("qLow"), " argument is not supported for CAGEset objects.")
  if (!is.null(qUp))
    stop(dQuote("qUp"), " argument is not supported for CAGEset objects.")
           
  gr <- CCdataframe2granges(object@consensusClusters)
  if (! is.null(sample))
    gr$tpm <- gr$score <- consensusClustersTpm(object)[,sample]
  gr
})

#' @rdname consensusClusters
#' @examples 
#' consensusClustersGR( exampleCAGEexp, sample = 2
#'                    , returnInterquantileWidth = TRUE
#'                    , qLow = 0.1, qUp = 0.9)

setMethod( "consensusClustersGR", "CAGEexp"
         , function (object, sample, returnInterquantileWidth, qLow, qUp) {
  if(is.null(experiments(object)$consensusClusters))
    stop("No consensus clusters found.  See ", sQuote("?aggregateTagClusters"), " on how to create them.")
  cc <- rowRanges(consensusClustersSE(object))
  if(!is.null(sample)) {
    if (!is.null(qLow))
      mcols(cc)[[paste0("q_", qLow)]] <-
        consensusClustersQuantile(object, sample, qLow)
    if (!is.null(qUp))
      mcols(cc)[[paste0("q_", qUp)]]  <-
        consensusClustersQuantile(object, sample, qUp)
    if (returnInterquantileWidth == TRUE) {
      if (is.null(qLow) | is.null(qUp))
        stop( "Set ", sQuote("qLow"), " and ", sQuote("qUp")
            , " to specify the quantile positions used to calculate width.")
      mcols(cc)[["interquantile_width"]] = mcols(cc)[[paste0("q_", qUp )]] -
                                           mcols(cc)[[paste0("q_", qLow)]] + 1
    }
    cc$tpm <- cc$score <- consensusClustersTpm(object)[,sample]
  }
  cc
})


#' @name consensusClustersSE
#' @rdname consensusClusters
#' @return `consensusClustersSE` returns the [`SummarizedExperiment`] stored
#' in the `consensusClusters` experiment slot of the CAGEexp object.
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
#' @title Get consensus clusters from CAGEr objects
#' 
#' @description Extracts the information on consensus clusters from a [`CAGEr`]
#'              object.
#' 
#' @param object A [`CAGEr`] object.
#' 
#' @param sample Optional. Label of the CAGE dataset (experiment, sample) for
#'        which to extract sample-specific information on consensus clusters.
#' 
#' @param returnInterquantileWidth Should the interquantile width of consensus
#'        clusters in specified sample be returned.  Used only when `sample`
#'        argument is specified, otherwise ignored.
#' 
#' @param qLow,qUp Position of which quantile should be used as a left (lower)
#'        or right (upper) boundary when calculating interquantile width.  Used
#'        only when `sample` argument is specified and
#'        `returnInterquantileWidth = TRUE`, otherwise ignored.
#' 
#' @return `consensusClusters` returns a `data.frame` with information on
#' consensus clusters, including genomic coordinates.  When `sample` argument is
#' NOT specified, total CAGE signal across all CAGE datasets (samples) is
#' returned in the `tpm` column.  When `sample` argument is specified, the `tpm`
#' column contains CAGE signal of consensus clusters in that specific sample.
#' When `returnInterquantileWidth = TRUE`, additional sample-specific information
#' is returned, including position of the dominant TSS, and interquantile width
#' of the consensus clusters in the specified sample.
#' 
#' @author Vanja Haberle
#' 
#' @seealso [`consensusClusters<-()`]
#' 
#' @family CAGEr accessor methods
#' @family CAGEr clusters functions
#' 
#' @examples
#' head(consensusClusters(exampleCAGEset))
#' head(consensusClusters(exampleCAGEset, sample = "sample2"))
#' 
#' cumulativeCTSSdistribution(exampleCAGEset, "consensusClusters") # Defaults in object do not fit
#' quantilePositions(exampleCAGEset, "consensusClusters")
#' head(consensusClusters(exampleCAGEset, sample = "sample2"
#'                       , returnInterquantileWidth = TRUE
#'                       , qLow = 0.1, qUp = 0.9))
#' 
#' head(consensusClusters(exampleCAGEexp))
#' consensusClustersGR(exampleCAGEexp, 2)
#' 
#' @export

setGeneric( "consensusClusters"
          , function( object, sample = NULL
                    , returnInterquantileWidth = FALSE, qLow = NULL, qUp = NULL) {
  validSamples(object, sample)
  standardGeneric("consensusClusters")
})

#' @rdname consensusClusters

setMethod( "consensusClusters", "CAGEr"
         , function (object, sample, returnInterquantileWidth, qLow, qUp) {
    
    # Get CCs specifically for each class
    if (inherits(object, "CAGEset")) {
      cc <- object@consensusClusters
    } else if (inherits(object, "CAGEexp")) {
      cc <- rowRanges(consensusClustersSE(object))
      cc$consensus.cluster <- names(cc)
      cc <- CCgranges2dataframe(cc)
    } else stop("Unsupported CAGEr class.")
    
    # Return raw without interquantiles if no sample requested
    if(is.null(sample))
      return(cc)
    
    # No support for per-sample output for CAGEexp objects 
    if(inherits(object, "CAGEexp"))
      stop( sQuote("sample"), " option not supported for CAGEexp objects.  Use "
          , sQuote("consensusClustersGR()"), " instead.")
           
    # Set score to expression levels in the selected sample
    cc$tpm <- consensusClustersTpm(object)[,sample]
    
    if(! returnInterquantileWidth) {
      # Return only clusters expressed in a given sample
      return(cc[cc$tpm > 0,])
    } else {
      # If interquantile width requested, check consistensy with other options.
      if (is.null(qLow) | is.null(qUp))
        stop("No quantiles specified! Please specify which quantile positions should be used to calculate width (qLow and qUp arguments)!")
      
      # Get dominant TSS position.
      cc.cumsum <- CTSScumulativesCC(object)[[sample]]
      cc$dominant_ctss <- sapply(cc.cumsum, .get.dominant.ctss, isCumulative = T)
      cc$dominant_ctss <- cc$start + cc$dominant_ctss
      
      # Return only clusters expressed in a given sample
      cc <- cc[cc$tpm > 0,]
      
      # Get dominant TSS expression score
      ctss <- CTSSnormalizedTpm(object)[,c("chr", "pos", "strand", sample)]
      pick <- match( paste(  cc$chr,   cc$dominant_ctss,   cc$strand)
                   , paste(ctss$chr, ctss$pos + 1,       ctss$strand))
      
      cc$tpm.dominant_ctss <- ctss[[sample]][pick]
      cc <- cc[,c("consensus.cluster", "chr", "start", "end", "strand", "dominant_ctss", "tpm", "tpm.dominant_ctss")]
      
      cc.w <- merge(consensusClustersQuantileLow(object, sample), consensusClustersQuantileUp(object, sample))
      cc.w <- cc.w[,c(1, which(colnames(cc.w) == paste0("q_", qLow)), which(colnames(cc.w) == paste0("q_", qUp)))]
      cc.w$interquantile_width <- cc.w[,3] - cc.w[,2] + 1
      cc <- merge(cc, cc.w, by.x = "consensus.cluster", by.y = "cluster", all.x = T)
    }
  return(cc)
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
#' @param sample A single sample name or number, or \code{NULL} (default) for all samples.
#' @param q A quantile.
NULL

#' @name consensusClustersQuantileLow
#' @rdname consensusClustersQuantile

setGeneric("consensusClustersQuantileLow", function(object, samples = NULL) {
  validSamples(object, samples)
  standardGeneric("consensusClustersQuantileLow")
})

#' @rdname consensusClustersQuantile

setMethod("consensusClustersQuantileLow", "CAGEset", function (object, samples) {
  result <- object@consensusClustersQuantileLow
  if (identical(result, list()))
    stop("Quantile positions not calculated yet! Run 'quantilePositions()' first !")
  if (is.null(samples)) result else result[[samples]]
})

#' @rdname consensusClustersQuantile

setMethod("consensusClustersQuantileLow", "CAGEexp", function (object, samples)
  stop( "Not supported for ", sQuote("CAGEexp"), " objects. "
      , "Use ", sQuote("consensusClustersQuantile()"), " instead."))


#' @name consensusClustersQuantileUp
#' @rdname consensusClustersQuantile

setGeneric("consensusClustersQuantileUp", function(object, samples = NULL) {
  validSamples(object, samples)
  standardGeneric("consensusClustersQuantileUp")
})

#' @rdname consensusClustersQuantile

setMethod("consensusClustersQuantileUp", "CAGEset", function (object, samples) {
  result <- object@consensusClustersQuantileUp
  if (identical(result, list()))
    stop("Quantile positions not calculated yet! Run 'quantilePositions()' first !")
  if (is.null(samples)) result else result[[samples]]
})

#' @rdname consensusClustersQuantile

setMethod("consensusClustersQuantileUp", "CAGEexp", function (object, samples)
  stop( "Not supported for ", sQuote("CAGEexp"), " objects. "
      , "Use ", sQuote("consensusClustersQuantile()"), " instead."))

#' @rdname consensusClustersQuantile

setGeneric("consensusClustersQuantile", function(object, sample = NULL, q) {
  validSamples(object, sample)
  standardGeneric("consensusClustersQuantile")
})

#' @rdname consensusClustersQuantile

setMethod("consensusClustersQuantile", "CAGEset", function (object, sample, q)
  stop( "Not supported for ", sQuote("CAGEexp"), " objects. "))

#' @rdname consensusClustersQuantile

setMethod("consensusClustersQuantile", "CAGEexp", function (object, sample, q) {
  if(is.null(assays(consensusClustersSE(object))[[paste0("q_", q)]]))
    stop("Quantile ", sQuote(q), " not found.")
  if (is.null(sample)) {
    assays(consensusClustersSE(object))[[paste0("q_", q)]]
  } else {
    assays(consensusClustersSE(object))[[paste0("q_", q)]][[sample]]
}})


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
#' @description Extracts a table with normalized CAGE tag values for consensus
#' clusters across all samples from a \code{\link{CAGEr}} object.
#' 
#' @param object A CAGEr object.
#' 
#' @return Returns the \code{matrix} of normalized expression values of CAGE clusters
#' across all samples.
#' 
#' @author Vanja Haberle
#' 
#' @family CAGEr clustering methods
#' @seealso \code{\link{consensusClusters}}
#' 
#' @examples
#' head(consensusClustersTpm(exampleCAGEset))
#' head(consensusClustersTpm(exampleCAGEexp))
#' 
#' @importFrom SummarizedExperiment assay
#' @export consensusClustersTpm

setGeneric("consensusClustersTpm", function(object) standardGeneric("consensusClustersTpm"))

#' @rdname consensusClustersTpm

setMethod("consensusClustersTpm", "CAGEset", function (object)
	object@consensusClustersTpmMatrix)

#' @rdname consensusClustersTpm

setMethod("consensusClustersTpm", "CAGEexp", function (object)
  assay(consensusClustersSE(object), "normalized"))


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

#' @name seqNameTotalsSE
#' 
#' @title Retreives the SummarizedExperiment containing chromosome expression totals.
#' 
#' @description Get or set a \code{SummarizedExperiment} summarising whole-chromosome
#' expression levels in the experiment slot \code{seqNameTotals} and the sample metadata
#' of the \code{\link{CAGEexp}} object.
#' 
#' @param object A \code{CAGEexp} object.
#' 
#' @family CAGEr accessor methods
#' @seealso summariseChrExpr
#' 
#' @author Charles Plessy
#' 
#' @examples 
#' seqNameTotalsSE(exampleCAGEexp)
#' 
#' @export

setGeneric("seqNameTotalsSE", function(object) standardGeneric("seqNameTotalsSE"))

#' @rdname seqNameTotalsSE

setMethod("seqNameTotalsSE", "CAGEset", function (object)
	stop("Not implemented for the CAGEset class."))

#' @rdname seqNameTotalsSE

setMethod("seqNameTotalsSE", "CAGEexp", function (object)
  experiments(object)$seqNameTotals)