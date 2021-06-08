#' @include CAGEr.R CTSS.R

####################################################################
# Functions for setting internal data in CAGEexp objects
#

#' @name genomeName<-
#' @rdname genomeName
#' @param value The name of a \code{BSgenome} package.
#' @family CAGEr setter methods
#' @author Charles Plessy
#' @export

setGeneric("genomeName<-", function(object, value) standardGeneric("genomeName<-"))

#' @rdname genomeName

setMethod("genomeName<-", "CAGEexp", function (object, value){
  metadata(object)$genomeName <- value
	if (validObject(object)) object
})

#' @rdname genomeName

setMethod("genomeName<-", "CTSS", function (object, value) {
  metadata(object)$genomeName <- value
	if (validObject(object)) object
})


#' @name inputFiles<-
#' @rdname inputFiles
#' @family CAGEr setter methods
#' @param value A character vector with one file path per sample.
#' @author Charles Plessy
#' @export

setGeneric("inputFiles<-", function(object, value) standardGeneric("inputFiles<-"))

#' @rdname inputFiles

setMethod("inputFiles<-", "CAGEexp", function (object, value){
  object$inputFiles <- value
  if (validObject(object)) object
})


#' @name inputFilesType<-
#' @rdname inputFilesType
#' @family CAGEr setter methods
#' @param value A character vector with one file type per sample.
#' @author Charles Plessy
#' @export

setGeneric("inputFilesType<-", function(object, value) standardGeneric("inputFilesType<-"))

#' @rdname inputFilesType

setMethod("inputFilesType<-", "CAGEexp", function (object, value){
  object$inputFilesType <- value
	if (validObject(object)) object
})


#' @name sampleLabels<-
#' @rdname sampleLabels
#' @family CAGEr setter methods
#' @param value A character vector with a unique and valid name for each sample.
#'        The \code{names} attributes indicate the colors.
#' @author Charles Plessy
#' @export

setGeneric("sampleLabels<-", function(object, value) standardGeneric("sampleLabels<-"))

#' @rdname sampleLabels

setMethod("sampleLabels<-", "CAGEexp", function (object, value){
  if (length(sampleLabels(object)) != length(value))
    stop("Number of labels differ from number of samples.")
  object$sampleLabels       <- unname(value)
  rownames(colData(object)) <- unname(value)
  if (! is.null(names(value)))
    object$Colors <- names(value)
  if (validObject(object)) object
})

#' @rdname sampleLabels

setMethod("sampleLabels<-", "CTSS", function (object, value){
	object@metadata$sampleLabels <- value
	if (validObject(object)) object
})

# librarySizes
#
# Not exported as it does not make sense to set library sizes after the data is loaded.

setGeneric("librarySizes<-", function(object, value) standardGeneric("librarySizes<-"))

setMethod("librarySizes<-", "CAGEexp", function (object, value){
  object$librarySizes <- value
  if (validObject(object)) object
})


#' @name CTSScoordinatesGR<-
#' @rdname CTSScoordinates
#' @param value Coordinates to update, in a format according to the function name.
#' @export

setGeneric("CTSScoordinatesGR<-", function(object, value) standardGeneric("CTSScoordinatesGR<-"))

#' @rdname CTSScoordinates

setMethod("CTSScoordinatesGR<-", "CAGEexp", function (object, value){
  if (! is(value, "GRanges")) stop("Value must be a GRanges object.")
  rowRanges(object@ExperimentList$tagCountMatrix) <- value
  if (validObject(object)) object
})


#' @name CTSStagCountSE<-
#' @rdname CTSScoordinates
#' @export

setGeneric("CTSStagCountSE<-", function(object, value) standardGeneric("CTSStagCountSE<-"))

#' @rdname CTSScoordinates

setMethod("CTSStagCountSE<-", "CAGEexp", function (object, value){
  if (! is(value, "RangedSummarizedExperiment"))
    stop("Value must be a RangedSummarizedExperiment object.")
  if (! all(colnames(value) == sampleLabels(object)))
    stop ("The CTSS data must match the CAGEexp object, with samples in the same order.")
  if (length(experiments(object)) == 0) {
    object <- MultiAssayExperiment( experiments = ExperimentList(tagCountMatrix=value)
                                  , colData = colData(object)
                                  , metadata = metadata(object))
    class(object) <- structure("CAGEexp", package = "CAGEr")
  } else if (is.null(object[["tagCountMatrix"]])) {
    object <- c(object, tagCountMatrix=value)
  } else {
    object[["tagCountMatrix"]] <- value
  }
  if (validObject(object)) object
})

#' @name filteredCTSSidx<-
#' 
#' @noRd 
#' 
#' @param value Logical

setGeneric( "filteredCTSSidx<-"
          , function(object, value) standardGeneric("filteredCTSSidx<-"))

setMethod("filteredCTSSidx<-", "CAGEexp", function (object, value) {
  CTSScoordinatesGR(object)$filteredCTSSidx <- value
  if (validObject(object)) object
})


#' @name CTSSclusteringMethod<-
#' 
#' @rdname CTSSclusteringMethod
#' 
#' @param value character

setGeneric( "CTSSclusteringMethod<-"
          , function(object, value) standardGeneric("CTSSclusteringMethod<-"))

#' @rdname CTSSclusteringMethod

setMethod("CTSSclusteringMethod<-", "GRangesList", function (object, value) {
  metadata(object)$clusteringMethod <- value
  if (validObject(object)) object
})

#' @rdname CTSSclusteringMethod

setMethod("CTSSclusteringMethod<-", "CAGEexp", function (object, value) {
  CTSSclusteringMethod(metadata(object)$tagClusters) <- value
  # extrat directly TCs from metadata slot because tagClustersGR does more that
  # is not needed here.
  if (validObject(object)) object
})


#' @name  CTSScumulativesTagClusters<-
#' 
#' @rdname CTSScumulativesTagClusters
#' 
#' @param value CTSScumulativesTagClusters data

setGeneric( "CTSScumulativesTagClusters<-"
          , function(object, value) standardGeneric("CTSScumulativesTagClusters<-"))

#' @rdname CTSScumulativesTagClusters

setMethod("CTSScumulativesTagClusters<-", "CAGEexp", function (object, value) {
  metadata(object)$CTSScumulativesTagClusters <- value
  if (validObject(object)) object
})

#' @name tagClustersGR<-
#' @rdname tagClusters
#' 
#' @param value A \code{\link{TagClusters}} object.

setGeneric( "tagClustersGR<-"
          , function(object, sample = NULL, value)
            standardGeneric("tagClustersGR<-"))

#' @rdname tagClusters

setMethod("tagClustersGR<-", c(object = "CAGEexp", value = "TagClusters"), function (object, sample, value) {
  validSamples(object, sample)
	metadata(object)$tagClusters[[sample]] <- value
  if (validObject(object)) object
})

#' @rdname tagClusters

setMethod("tagClustersGR<-", c("CAGEexp", "missing", "GRangesList"), function (object, sample, value) {
	metadata(object)$tagClusters <- value
  if (validObject(object)) object
})

#' @name consensusClusters<-
#' @rdname consensusClusters-set
#' 
#' @title Set consensus clusters from CAGEr objects
#' 
#' @description Set the information on consensus clusters in a [`CAGEr`]
#'              object.
#' 
#' @param object A [`CAGEr`] object.
#' @param value A \code{data.frame} of consensus clusters
#' 
#' @details These setter methods are mostly for internal use, but are exported
#' in case they may be useful to advanced users.
#' 
#' @author Vanja Haberle
#' @author Charles Plessy
#' 
#' @aliases consensusClustersSE<-
#' @export

setGeneric("consensusClustersSE<-", function(object, value) standardGeneric("consensusClustersSE<-"))

#' @rdname consensusClusters-set

setMethod( "consensusClustersSE<-"
         , c("CAGEexp", "RangedSummarizedExperiment")
         , function (object, value) {
  if (! all(colnames(value) == sampleLabels(object)))
    stop ("The expression data must match the CAGEexp object, with samples in the same order.")
  if (is.null(object[["consensusClusters"]])) {
    object <- c(object, consensusClusters = value)
  } else {
    object[["consensusClusters"]] <- value
  }
  if (validObject(object)) object
})

#' @name consensusClustersGR<-
#' @rdname consensusClusters-set
#' @export

setGeneric("consensusClustersGR<-", function(object, value) standardGeneric("consensusClustersGR<-"))

#' @rdname consensusClusters-set

setMethod("consensusClustersGR<-", "CAGEexp", function (object, value){
  if (! is(value, "GRanges")) stop("Value must be a GRanges object.")
  rowRanges(object@ExperimentList$consensusClusters) <- value
  if (validObject(object)) object
})


#' @name consensusClustersQuantileLow<-
#' @rdname consensusClustersQuantile

setGeneric( "consensusClustersQuantileLow<-"
          , function(object, samples = NULL, value)
              standardGeneric("consensusClustersQuantileLow<-"))

#' @name consensusClustersQuantileUp<-
#' @rdname consensusClustersQuantile

setGeneric( "consensusClustersQuantileUp<-"
          , function(object, samples = NULL, value)
              standardGeneric("consensusClustersQuantileUp<-"))

#' @name `CTSScumulativesCC<-`
#' @noRd

setGeneric("CTSScumulativesCC<-", function(object, value) standardGeneric("CTSScumulativesCC<-"))

setMethod("CTSScumulativesCC<-", "CAGEexp", function (object, value){
	metadata(object)$CTSScumulativesConsensusClusters <- value
  if (validObject(object)) object
})


# GeneExpSE
# 
# Since the SummarizedExperiment can hold normalized and non-normalized values,
# let's name it "GeneExp" instead of "GeneTagCount" if we would follow the
# historical CTSS name pattern of CAGEset objects.
#
# Not exported for the moment.

setGeneric("GeneExpSE<-", function(object, value) standardGeneric("GeneExpSE<-"))

setMethod("GeneExpSE<-", "CAGEexp", function (object, value){
  if (! is(value, "SummarizedExperiment"))
    stop("Value must be a SummarizedExperiment object.")
  if (is(value, "RangedSummarizedExperiment"))
    stop("Value must not be a RangedSummarizedExperiment object. ",
         "(Gene symbols have no ranged coordinates).")
  if (! all(colnames(value) == sampleLabels(object)))
    stop ("The CTSS data must match the CAGEexp object, with samples in the same order.")
  if (is.null(object[["geneExpMatrix"]])) {
    object <- c(object, geneExpMatrix = value)
  } else {
    object[["geneExpMatrix"]] <- value
  }
  if (validObject(object)) object
})

#' @rdname seqNameTotalsSE
#' @param value A SummarizedExperiment object where rows represent reference sequences
#'        such as chromosomes.

setGeneric("seqNameTotalsSE<-", function(object, value) standardGeneric("seqNameTotalsSE<-"))

setMethod( "seqNameTotalsSE<-"
         , c("CAGEexp", "SummarizedExperiment")
         , function (object, value) {
  if (! all(colnames(value) == sampleLabels(object)))
    stop ("The expression data must match the CAGEexp object, with samples in the same order.")
  if (is.null(object[["seqNameTotals"]])) {
    object <- c(object, seqNameTotals = value)
  } else {
    object[["seqNameTotals"]] <- value
  }
  if (validObject(object)) object
})

#' @name setColors
#' 
#' @title Set colors for samples
#' 
#' @description Assigns one color to each sample in the CAGEr object.  These
#' colors are used in various plots and exported tracks to consistently
#' represent corresponding samples.
#' 
#' @param object A \code{\link{CAGEr}} object.
#' @param colors A character vector of one valid \R color specification per
#'   sample (see \code{\link{col2rgb}} for details).  Provided colors are
#'   assigned to samples in the order they are returned by the
#'   \code{\link{sampleLabels}} function.
#' 
#' @return Assigns one color to each sample in the CAGEr object and modifies it
#' in place.
#' 
#' @author Vanja Haberle
#' 
#' @family CAGEr setter methods
#' 
#' @importFrom grDevices col2rgb
#' @importFrom grDevices rgb
#' 
#' @examples
#' 
#' sampleLabels(exampleCAGEexp)
#' setColors(exampleCAGEexp, 5)
#' sampleLabels(exampleCAGEexp)
#' setColors(exampleCAGEexp, c("#ff0000ff", "#CCFF00", "blue", "grey", 1))
#' sampleLabels(exampleCAGEexp)
#' setColors(exampleCAGEexp, c("red", "darkgreen", "blue", "grey", "black"))
#' sampleLabels(exampleCAGEexp)
#' 
#' @export

setGeneric("setColors", function(object, colors = NULL) standardGeneric("setColors"))

#' @rdname setColors

setMethod("setColors", "CAGEr", function (object, colors){

  sample.labels <- sampleLabels(object)
  
  if(length(colors) == 1 & is.numeric(colors)){
  	names(sample.labels) <- rainbow(n = length(sample.labels))
  }else if(length(colors) != length(sample.labels)){
  	stop(paste("Number of provided colors must match the number of samples in the CAGEr object, i.e. must be ", length(sample.labels), "!", sep = ""))
  }else{
    names(sample.labels) <- sapply(colors, function(x){
      rgb.col <- tryCatch( col2rgb(x, alpha = TRUE)
                         , error = function(e) stop(dQuote(x), " is not a valid color. See col2rgb() for details.", call. = FALSE))
      do.call(rgb, c(as.list(rgb.col), maxColorValue = 255))
    })
  }
  
  sampleLabels(object) <- sample.labels
  object
})
