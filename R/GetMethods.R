#' @include AggregationFunctions.R ClusteringFunctions.R CAGEr.R CTSS.R

################################################################
# Functions for retrieving data from CAGEexp objects

#' @name genomeName
#' 
#' @title Extracting genome name from CAGEr objects
#' 
#' @description Extracts the name of a referent genome from a
#' `CAGEexp` or a `CTSS` object.
#' 
#' @param object A CAGEexp or a CTSS object.
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
#' genomeName(exampleCAGEexp)
#' 
#' @export

setGeneric("genomeName", function(object) standardGeneric("genomeName"))

#' @rdname genomeName

setMethod("genomeName", "CAGEexp", function (object) metadata(object)$genomeName)

#' @rdname genomeName

setMethod("genomeName", "CTSS", function (object) metadata(object)$genomeName)


#' @name inputFiles
#' 
#' @title Extracting paths to input files from CAGEr objects
#' 
#' @description Extracts the paths to CAGE data input files from
#' code{\link{CAGEexp}} objects.
#' 
#' @param object A CAGEexp object.
#' 
#' @return Returns a character vector of paths to CAGE data input files.
#' 
#' @family CAGEr accessor methods
#' 
#' @author Vanja Haberle
#' 
#' @examples 
#' inputFiles(exampleCAGEexp)
#' 
#' @export

setGeneric(
name="inputFiles",
def=function(object){
	standardGeneric("inputFiles")
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
#' files from \code{\link{CAGEexp}} objects.
#' 
#' @param object A CAGEexp object.
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
#' inputFilesType(exampleCAGEexp)
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
signature(object = "CAGEexp"),
function (object){
  object$inputFilesType
})


#' @name librarySizes
#' 
#' @title Extracting library sizes from CAGEr objects
#' 
#' @description Extracts the library sizes (total number of CAGE tags) for all CAGE datasets
#' from \code{\link{CAGEexp}} objects.
#' 
#' @param object A CAGEexp object.
#' 
#' @details Library sizes are calculated when loading data with the \code{getCTSS}
#' function and stored in the \code{librarySizes} column of the \code{colData} of
#' \code{CAGEexp} objects.
#' 
#' @return Returns an integer vector of total number of CAGE tags (library size) for all CAGE
#' datasets in the CAGEr object.
#' 
#' @seealso \code{\link{getCTSS}}
#' 
#' @examples 
#' librarySizes(exampleCAGEexp)
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
signature(object = "CAGEexp"),
function (object){
  as.integer(object$librarySizes)
})

#' @rdname CTSScoordinates
#' @name CTSScoordinatesGR
#' 
#' @title Genomic coordinates of TSSs from a `CAGEr` object
#' 
#' @description Extracts the genomic coordinates of all detected TSSs
#' from [CAGEexp] objects.
#' 
#' @param object A `CAGEexp` object.
#' 
#' @seealso
#' [`getCTSS`]
#' 
#' @examples
#' CTSScoordinatesGR(exampleCAGEexp)
#'
#' @author Vanja Haberle
#' @author Charles Plessy
#' @family CAGEr accessor methods
#' @export
#' 
#' @return `CTSScoordinatesGR` returns the coordinates as a [CTSS()] object
#' wrapping genomic ranges.  A `filteredCTSSidx` column metadata will be present
#' if [clusterCTSS()] was ran earlier.
#' 
#' @seealso [`clusterCTSS`]
#' 
#' 
#' @importFrom GenomeInfoDb genome genome<-
#' @importFrom IRanges IRanges
#' 
#' @examples
#' CTSScoordinatesGR(exampleCAGEexp)
#' 
#' @export

setGeneric("CTSScoordinatesGR", function(object) standardGeneric("CTSScoordinatesGR"))

#' @rdname CTSScoordinates

setMethod("CTSScoordinatesGR", "CAGEexp", function (object)
  rowRanges(CTSStagCountSE(object)))

#' @rdname CTSStagCount
#' @name CTSStagCountDF
#' 
#' @title Raw CAGE TSSs expression counts
#' 
#' @description Extracts the tag count for all detected TSSs in all CAGE datasets
#'              from [`CAGEexp`] objects.
#' 
#' @param object A `CAGEexp` object.
#' @param samples For `CTSStagCountGR` only: name(s) or number(s) identifying
#' sample(s) or "all" to return a `GRangesList` of all the samples.
#'  
#' @return Returns an object with number of CAGE tags supporting each TSS
#' (rows) in every CAGE dataset (columns).  The class of the object depends on the
#' function being called:
#' 
#' * `CTSStagCountDF`: A [`DataFrame`] of [`Rle`] integers.
#' * `CTSStagCountDA`: A [`DelayedArray`] wrapping a `DataFrame` of `Rle` integers.
#' * `CTSStagCountSE`: A [`RangedSummarizedExperiment`]` containing a `DataFrame`
#'    of `Rle` integers.
#' * `CTSStagCountGR`: A `CTSS` object (wrapping `GRanges`) containing a `score`
#'    column indicating expression values for a given sample, or a
#'   `GRangesList` of `CTSS` objects.
#' 
#' @seealso [getCTSS()]
#' 
#' @author Vanja Haberle
#' @author Charles Plessy
#' 
#' @family CAGEr accessor methods
#' @family CAGEr CTSS methods
#' 
#' @examples 
#' CTSStagCountDF(exampleCAGEexp)
#'  
#' @export

setGeneric("CTSStagCountDF", function(object) standardGeneric("CTSStagCountDF"))

#' @rdname CTSStagCount

setMethod("CTSStagCountDF", "CAGEexp", function (object)
  assay(CTSStagCountSE(object)))


#' @name CTSStagCountDA
#' @rdname CTSStagCount
#' 
#' @examples 
#' CTSStagCountDA(exampleCAGEexp)
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
#' @examples 
#' CTSStagCountGR(exampleCAGEexp, 1)
#' CTSStagCountGR(exampleCAGEexp, "all")
#'  
#' @export

setGeneric("CTSStagCountGR", function(object, samples) {
  validSamples(object, samples)
  standardGeneric("CTSStagCountGR")
})

#' @rdname CTSStagCount

setMethod( "CTSStagCountGR", "CAGEexp", function (object, samples) {
  if (samples == "all") {
    l <- lapply(seq_along(sampleLabels(object)), CTSStagCountGR, object = object)
    return(GRangesList(l))
  }
  if (is.character(samples)) samples <- which(sampleLabels(object) == samples)
  gr <- CTSScoordinatesGR(object)
  score(gr) <- CTSStagCountDF(object)[[samples]]
  gr <- gr[score(gr) != 0]
  sampleLabels(gr) <- sampleLabels(object)[samples]
  gr
})

#' @name CTSStagCountSE
#' @rdname CTSStagCount
#' 
#' @examples 
#' CTSStagCountSE(exampleCAGEexp)
#' 
#' @export

setGeneric("CTSStagCountSE", function(object) standardGeneric("CTSStagCountSE"))

#' @rdname CTSStagCount

setMethod("CTSStagCountSE", "CAGEexp", function (object) {
  se <- experiments(object)$tagCountMatrix
  if (is.null(se)) stop("Could not find CTSS tag counts, see ", sQuote("?getCTSS"), ".")
  se
})

#' @name CTSSnormalizedTpmDF
#' @rdname CTSSnormalizedTpm
#' 
#' @title Extracting normalized CAGE signal for TSSs from CAGEr objects
#' 
#' @description Extracts the normalized CAGE signal for all detected TSSs
#' in all CAGE datasets from [`CAGEexp`] objects.
#' 
#' @param object A `CAGEexp` object.
#' 
#' @seealso \code{\link{normalizeTagCount}}
#' 
#' @examples 
#' CTSSnormalizedTpmDF(exampleCAGEexp)
#' 
#' @author Vanja Haberle
#' @author Charles Plessy
#' @family CAGEr accessor methods
#' @export
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
#' @examples 
#' CTSSnormalizedTpmGR(exampleCAGEexp, 1)
#' exampleCAGEexp |> CTSSnormalizedTpmGR("all") 
#'  
#' @export

setGeneric("CTSSnormalizedTpmGR", function(object, samples) {
  validSamples(object, samples)
  standardGeneric("CTSSnormalizedTpmGR")
  })

#' @rdname CTSSnormalizedTpm

setMethod( "CTSSnormalizedTpmGR", "CAGEexp", function (object, samples) {
  if (samples == "all") {
    l <- lapply(seq_along(sampleLabels(object)), CTSSnormalizedTpmGR, object = object)
    return(GRangesList(l))
  }
  if (is.character(samples)) samples <- which(sampleLabels(object) == samples)
  gr <- CTSScoordinatesGR(object)
  score(gr) <- CTSSnormalizedTpmDF(object)[[samples]]
  gr <- gr[score(gr) != 0]
  sampleLabels(gr) <- sampleLabels(object)[samples]
  gr
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
#' CTSSclusteringMethod(exampleCAGEexp)
#' 
#' @export CTSSclusteringMethod

setGeneric("CTSSclusteringMethod", function(object) standardGeneric("CTSSclusteringMethod"))

#' @rdname CTSSclusteringMethod

setMethod("CTSSclusteringMethod", "GRangesList", function (object)
	metadata(object)$clusteringMethod)

#' @rdname CTSSclusteringMethod

setMethod("CTSSclusteringMethod", "CAGEexp", function (object)
  CTSSclusteringMethod(metadata(object)$tagClusters))


#' @name tagClustersGR
#' @rdname tagClusters
#' 
#' @title Extract tag clusters (TCs) for individual CAGE experiments
#' 
#' @description Extracts tag clusters (TCs) produced by [`clusterCTSS`] function
#' for a specified CAGE experiment from a [`CAGEexp`] object.
#' 
#' @param object A `CAGEexp` object.
#' 
#' @param sample Label of the CAGE dataset (experiment, sample) for which to
#' extract tag clusters. If `samples = NULL`, a list of all the clusters for
#' each sample is returned.
#' 
#' @param returnInterquantileWidth Return the interquantile width for each tag cluster.
#' 
#' @param qLow,qUp Position of which quantile should be used as a left (lower)
#' or right (upper) boundary (for `qLow` and `qUp` respectively) when
#' calculating interquantile width.  Default value `NULL` results in using the
#' start coordinate of the cluster.  Used only when
#' `returnInterquantileWidth = TRUE`, otherwise ignored.
#' 
#' @return Returns a `GRangesList` or a `GRanges` object with genomic coordinates,
#' position of dominant TSS, total CAGE signal and additional information for
#' all TCs from specified CAGE dataset (sample).  If
#' `returnInterquantileWidth = TRUE`, interquantile width for each TC is also
#' calculated using provided quantile positions.
#' 
#' @author Vanja Haberle
#' @author Charles Plessy
#' 
#' @family CAGEr accessor methods
#' @family CAGEr clusters functions
#' @export
#' 
#' @examples
#' tagClustersGR( exampleCAGEexp, "Zf.high", TRUE, 0.1, 0.9 )
#' tagClustersGR( exampleCAGEexp, 1
#'              , returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9 )
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

setMethod("filteredCTSSidx", "CAGEexp", function (object){
  rowData(CTSStagCountSE(object))$filteredCTSSidx
})


#' @name consensusClustersGR
#' @rdname consensusClusters
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
#' @return `consensusClustersGR` returns a [`ConsensusClusters`] object, which
#' wraps the [`GRanges`] class.  The `score` columns indicates the
#' normalised expression value of each cluster, either across all samples
#' (`sample = NULL`), or for the selected sample.  The legacy `tpm` column may
#' be removed in the future.  When `sample` argument is
#' NOT specified, total CAGE signal across all CAGE datasets (samples) is
#' returned in the `tpm` column.  When `sample` argument is specified, the `tpm`
#' column contains CAGE signal of consensus clusters in that specific sample.
#' When `returnInterquantileWidth = TRUE`, additional sample-specific information
#' is returned, including position of the dominant TSS, and interquantile width
#' of the consensus clusters in the specified sample.
#' 
#' @author Vanja Haberle
#' @author Charles Plessy
#' 
#' @seealso [`consensusClusters<-()`]
#' 
#' @family CAGEr accessor methods
#' @family CAGEr clusters functions
#' 
#' @examples
#' consensusClustersGR( exampleCAGEexp, sample = 2
#'                    , returnInterquantileWidth = TRUE
#'                    , qLow = 0.1, qUp = 0.9)
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

setMethod( "consensusClustersGR", "CAGEexp"
         , function (object, sample, returnInterquantileWidth, qLow, qUp) {
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

setMethod("consensusClustersSE", "CAGEexp", function (object) {
  if(is.null(experiments(object)$consensusClusters))
    stop("No consensus clusters found.  See ", sQuote("?aggregateTagClusters"), " on how to create them.")
  experiments(object)$consensusClusters
})


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

setMethod( "consensusClustersDESeq2", "CAGEexp"
         , function (object, design) {
  if (! requireNamespace("DESeq2"))
    stop("This function requires the ", dQuote("DESeq2"), " package; please install it.")
  DESeq2::DESeqDataSetFromMatrix( countData = assay(consensusClustersSE(object))
                                , colData   = colData(object)
                                , rowData   = consensusClustersGR(object)
                                , design    = design)
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

setMethod("consensusClustersQuantileUp", "CAGEexp", function (object, samples)
  stop( "Not supported for ", sQuote("CAGEexp"), " objects. "
      , "Use ", sQuote("consensusClustersQuantile()"), " instead."))

#' @rdname consensusClustersQuantile

setGeneric("consensusClustersQuantile", function(object, sample = NULL, q) {
  validSamples(object, sample)
  standardGeneric("consensusClustersQuantile")
})

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
#' @param object A [`CAGEexp`] object.
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
#' @title Extracting consensus clusters tpm matrix from CAGEr object
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
#' @seealso \code{\link{consensusClustersSE}}
#' 
#' @examples
#' head(consensusClustersTpm(exampleCAGEexp))
#' 
#' @importFrom SummarizedExperiment assay
#' @export consensusClustersTpm

setGeneric("consensusClustersTpm", function(object) standardGeneric("consensusClustersTpm"))

#' @rdname consensusClustersTpm

setMethod("consensusClustersTpm", "CAGEexp", function (object)
  assay(consensusClustersSE(object), "normalized"))


#' @name expressionClasses
#' 
#' @title Extract labels of _expression classes_
#' 
#' @description Retrieves labels of _expression classes_ of individual CTSSs
#' or consensus clusters from a `CAGEr` object.
#' 
#' @param object A [`CAGEr`] object.
#' 
#' @return Returns a [`Rle`]-encoded vector of labels of _expression classes_.
#' The number of labels matches the number of expression clusters returned by
#' [`getExpressionProfiles`] function.
#' 
#' @family CAGEr expression clustering functions
#' @family CAGEr accessor methods
#' 
#' @examples
#' expressionClasses(CTSScoordinatesGR(exampleCAGEexp))
#' exampleCAGEexp |> consensusClustersGR() |> expressionClasses()
#' 
#' @export

setGeneric( "expressionClasses", function(object)
  standardGeneric("expressionClasses"))

#' @rdname expressionClasses

setMethod("expressionClasses", "CTSS", function (object) {
  classes <- object$exprClass
    if (is.null(classes)) stop("No expression clustering of CTSSs has been done yet!")
  classes
})

#' @rdname expressionClasses

setMethod("expressionClasses", "ConsensusClusters", function (object) {
  classes <- object$exprClass
  if (is.null(classes)) stop("No expression clustering of consensus clusters has been done yet!")
  classes
})


#' @name CTSSexpressionClusteringMethod
#' @noRd 

setGeneric("CTSSexpressionClusteringMethod", function(object) 
  standardGeneric("CTSSexpressionClusteringMethod"))

setMethod("CTSSexpressionClusteringMethod", "CAGEexp", function (object)
  metadata(object)$CTSSexpressionClusteringMethod)


#' @name consensusClustersExpressionClusteringMethod
#' @noRd 

setGeneric("consensusClustersExpressionClusteringMethod", function(object) 
  standardGeneric("consensusClustersExpressionClusteringMethod"))

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

setMethod("seqNameTotalsSE", "CAGEexp", function (object)
  experiments(object)$seqNameTotals)
