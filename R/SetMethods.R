#' @include CAGEr.R

####################################################################
# Functions for setting internal data in CAGEset and CAGEexp objects
#

#' @name genomeName<-
#' @rdname genomeName
#' @param value The name of a \code{BSgenome} package.
#' @family CAGEr setter methods
#' @author Charles Plessy
#' @export

setGeneric("genomeName<-", function(object, value) standardGeneric("genomeName<-"))

#' @rdname genomeName

setMethod("genomeName<-", "CAGEset", function (object, value){
	object@genomeName <- value
	if (validObject(object)) object
})

#' @rdname genomeName

setMethod("genomeName<-", "CAGEexp", function (object, value){
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

setMethod("inputFiles<-", "CAGEset", function (object, value){
	object@inputFiles <- value
	if (validObject(object)) object
})

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

setMethod("inputFilesType<-", "CAGEset", function (object, value){
	object@inputFilesType <- value
	if (validObject(object)) object
})

#' @rdname inputFilesType

setMethod("inputFilesType<-", "CAGEexp", function (object, value){
  object$inputFilesType <- value
	if (validObject(object)) object
})


#' @name sampleLabels<-
#' @rdname sampleLabels
#' @family CAGEr setter methods
#' @param value A character vector with a unique and valid name for each sample.
#'        the \code{names} attributes indicate the colors.
#' @author Charles Plessy
#' @export

setGeneric("sampleLabels<-", function(object, value) standardGeneric("sampleLabels<-"))

#' @rdname sampleLabels

setMethod("sampleLabels<-", "CAGEset", function (object, value){
	object@sampleLabels <- value
	if (validObject(object)) object
})

#' @rdname sampleLabels

setMethod("sampleLabels<-", "CAGEexp", function (object, value){
  if (length(sampleLabels(object)) != length(value))
    stop("Number of labels differ from number of samples.")
  object$sampleLabels       <- value
  rownames(colData(object)) <- value
  if (validObject(object)) object
})


# librarySizes
#
# Not exported as it does not make sense to set library sizes after the data is loaded.

setGeneric("librarySizes<-", function(object, value) standardGeneric("librarySizes<-"))

setMethod("librarySizes<-", "CAGEset", function (object, value){
	object@librarySizes <- value
	if (validObject(object)) object
})

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

setMethod("CTSScoordinatesGR<-", "CAGEset", function (object, value){
	stop("Not implemented for the CAGEset class.")
})

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

setMethod("CTSStagCountSE<-", "CAGEset", function (object, value){
	stop("Not implemented for the CAGEset class.")
})

#' @rdname CTSScoordinates

setMethod("CTSStagCountSE<-", "CAGEexp", function (object, value){
  if (! is(value, "RangedSummarizedExperiment"))
    stop("Value must be a RangedSummarizedExperiment object.")
  if (! all(colnames(value) == sampleLabels(object)))
    stop ("The CTSS data must match the CAGEexp object, with samples in the same order.")
  sampleMapSE <-
    listToMap(list(tagCountMatrix = data.frame( primary = sampleLabels(object)
                                              , colname = colnames(value))))
  sampleMap(object) <-
    rbind( sampleMap(object)[sampleMap(object)$assay != "tagCountMatrix",]
         , sampleMapSE)
  experiments(object)$tagCountMatrix <- value
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

setMethod("CTSScumulativesTagClusters<-", "CAGEset", function (object, value) {
	object@CTSScumulativesTagClusters <- value
	if (validObject(object)) object
})

#' @rdname CTSScumulativesTagClusters

setMethod("CTSScumulativesTagClusters<-", "CAGEexp", function (object, value) {
  metadata(object)$CTSScumulativesTagClusters <- value
  if (validObject(object)) object
})


#' @name tagClustersQuantileLow<-
#' @rdname tagClustersQuantile

setGeneric( "tagClustersQuantileLow<-"
          , function(object, value, samples = NULL)
              standardGeneric("tagClustersQuantileLow<-"))

#' @rdname tagClustersQuantile

setMethod("tagClustersQuantileLow<-", "CAGEset", function (object, value, samples) {
  validSamples(object, samples)
  if(is.null(samples)) {
	  object@tagClustersQuantileLow <- value
  } else {
    object@tagClustersQuantileLow[[samples]] <- value
  }
	if (validObject(object)) object
})

#' @rdname tagClustersQuantile

setMethod("tagClustersQuantileLow<-", "CAGEexp", function (object, value, samples)
  stop("Not supported for CAGEexp.  Use tagClustersGR<- instead."))

#' @name tagClustersQuantileUp<-
#' @rdname tagClustersQuantile
#' 
setGeneric( "tagClustersQuantileUp<-"
          , function(object, value, samples = NULL)
              standardGeneric("tagClustersQuantileUp<-"))

#' @rdname tagClustersQuantile

setMethod("tagClustersQuantileUp<-", "CAGEset", function (object, value, samples) {
  validSamples(object, samples)
  if(is.null(samples)) {
	  object@tagClustersQuantileUp <- value
  } else {
    object@tagClustersQuantileUp[[samples]] <- value
  }
	if (validObject(object)) object
})

#' @rdname tagClustersQuantile

setMethod("tagClustersQuantileLow<-", "CAGEexp", function (object, value, samples)
  stop("Not supported for CAGEexp.  Use tagClustersGR<- instead."))


#' @name tagClustersGR<-
#' @rdname tagClusters
#' 
#' @param value A \code{\link{TagClusters}} object.

setGeneric( "tagClustersGR<-"
          , function(object, value, samples = NULL)
            standardGeneric("tagClustersGR<-"))

#' @rdname tagClusters

setMethod("tagClustersGR<-", "CAGEset", function (object, value, samples)
  stop("Not supported for CAGEset. Set tagClustersQuantileLow/Up instead."))

#' @rdname tagClusters

setMethod("tagClustersGR<-", c("CAGEexp", "TagClusters"), function (object, value, samples) {
  validSamples(object, samples)
	metadata(object)$tagClusters[[samples]] <- value
  if (validObject(object)) object
})

setMethod("tagClustersGR<-", c("CAGEexp", "GRangesList"), function (object, value, samples) {
	metadata(object)$tagClusters <- value
  if (validObject(object)) object
})

#' @name consensusClusters<-
#' @rdname consensusClusters
#' 
#' @param value A \code{data.frame} of consensus clusters

setGeneric("consensusClusters<-", function(object, value) standardGeneric("consensusClusters<-"))

#' @rdname consensusClusters

setMethod("consensusClusters<-", "CAGEset", function (object, value){
	object@consensusClusters <- value
  if (validObject(object)) object
})

#' @rdname consensusClusters

setMethod("consensusClusters<-", "CAGEexp", function (object, value){
  stop("Not supported for CAGEexp.")
})

#' @name consensusClustersSE<-
#' @rdname consensusClusters
#' @export

setGeneric("consensusClustersSE<-", function(object, value) standardGeneric("consensusClustersSE<-"))

#' @rdname consensusClusters

setMethod("consensusClustersSE<-", "CAGEset", function (object, value){
	stop("Not implemented for the CAGEset class.")
})

#' @rdname consensusClusters

setMethod("consensusClustersSE<-", c("CAGEexp", "RangedSummarizedExperiment"), function (object, value){
  if (! all(colnames(value) == sampleLabels(object)))
    stop ("The CTSS data must match the CAGEexp object, with samples in the same order.")
  sampleMapSE <-
    listToMap(list(consensusClusters = data.frame( primary = sampleLabels(object)
                                              , colname = colnames(value))))
  sampleMap(object) <-
    rbind( sampleMap(object)[sampleMap(object)$assay != "consensusClusters",]
         , sampleMapSE)
  experiments(object)$consensusClusters <- value
  if (validObject(object)) object
})

#' @name consensusClustersQuantileLow<-
#' @rdname consensusClustersQuantile

setGeneric("consensusClustersQuantileLow<-", function(object, value) standardGeneric("consensusClustersQuantileLow<-"))

#' @rdname consensusClustersQuantile

setMethod("consensusClustersQuantileLow<-", "CAGEset", function (object, value){
	object@consensusClustersQuantileLow <- value
	if (validObject(object)) object
})

#' @rdname consensusClustersQuantile

setMethod("consensusClustersQuantileLow<-", "CAGEexp", function (object, value){
  metadata(object)$consensusClustersQuantileLow <- value
  if (validObject(object)) object
})


#' @name consensusClustersQuantileUp<-
#' @rdname consensusClustersQuantile

setGeneric("consensusClustersQuantileUp<-", function(object, value) standardGeneric("consensusClustersQuantileUp<-"))

#' @rdname consensusClustersQuantile

setMethod("consensusClustersQuantileUp<-", "CAGEset", function (object, value){
	object@consensusClustersQuantileUp <- value
  if (validObject(object)) object
})

#' @rdname consensusClustersQuantile

setMethod("consensusClustersQuantileUp<-", "CAGEexp", function (object, value){
  metadata(object)$consensusClustersQuantileUp <- value
  if (validObject(object)) object
})



#' @name `CTSScumulativesCC<-`
#' @noRd

setGeneric("CTSScumulativesCC<-", function(object, value) standardGeneric("CTSScumulativesCC<-"))

setMethod("CTSScumulativesCC<-", "CAGEset", function (object, value){
	object@CTSScumulativesConsensusClusters <- value
  if (validObject(object)) object
})

setMethod("CTSScumulativesCC<-", "CAGEexp", function (object, value){
	metadata(object)$CTSScumulativesConsensusClusters <- value
  if (validObject(object)) object
})

#' @name consensusClustersGR<-
#' @rdname consensusClusters
#' @export

setGeneric("consensusClustersGR<-", function(object, value) standardGeneric("consensusClustersGR<-"))

#' @rdname consensusClusters

setMethod("consensusClustersGR<-", "CAGEset", function (object, value){
	stop("Not implemented for the CAGEset class.")
})

#' @rdname consensusClusters

setMethod("consensusClustersGR<-", "CAGEexp", function (object, value){
  if (! is(value, "GRanges")) stop("Value must be a GRanges object.")
  rowRanges(object@ExperimentList$consensusClusters) <- value
  if (validObject(object)) object
})


#' @name consensusClustersSE<-
#' @rdname consensusClusters

setGeneric("consensusClustersSE<-", function(object, value) standardGeneric("consensusClustersSE<-"))

#' @rdname consensusClusters

setMethod("consensusClustersSE<-", "CAGEset", function (object, value){
	stop("Not implemented for the CAGEset class.")
})

#' @rdname consensusClusters

setMethod( "consensusClustersSE<-"
         , c("CAGEexp", "RangedSummarizedExperiment")
         , function (object, value) {
  if (! all(colnames(value) == sampleLabels(object)))
    stop ("The expression data must match the CAGEexp object, with samples in the same order.")
  sampleMapSE <-
    listToMap(list(consensusClusters = data.frame( primary = sampleLabels(object)
                                                 , colname = colnames(value))))
  sampleMap(object) <-
    rbind( sampleMap(object)[sampleMap(object)$assay != "consensusClusters",]
         , sampleMapSE)
  experiments(object)$consensusClusters <- value
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

setMethod("GeneExpSE<-", "CAGEset", function (object, value){
	stop("Not implemented for the CAGEset class.")
})

setMethod("GeneExpSE<-", "CAGEexp", function (object, value){
  if (! is(value, "SummarizedExperiment"))
    stop("Value must be a SummarizedExperiment object.")
  if (is(value, "RangedSummarizedExperiment"))
    stop("Value must not be a RangedSummarizedExperiment object. ",
         "(Gene symbols have no ranged coordinates).")
  if (! all(colnames(value) == sampleLabels(object)))
    stop ("The CTSS data must match the CAGEexp object, with samples in the same order.")
  sampleMapSE <-
    listToMap(list(geneExpMatrix = data.frame( primary = sampleLabels(object)
                                             , colname = colnames(value))))
  sampleMap(object) <-
    rbind( sampleMap(object)[sampleMap(object)$assay != "geneExpMatrix",]
         , sampleMapSE)
  experiments(object)$geneExpMatrix <- value
  if (validObject(object)) object
})

#' @name setColors
#' 
#' @title Setting colors for samples
#' 
#' @description Assigns one color to each sample in the CAGEr object.  These colors are
#' used in various plots and exported tracks to consistently represent corresponding samples.
#' 
#' @param object A \code{\link{CAGEr}} object.
#' @param colors A character vector of valid color names.  For a complete list of valid color
#'               names see the \code{\link{colors}} function.  Alternatively, it can be a
#'               character vector of colors specified in hexadecimal format (\emph{e.g.}
#'               "#FF0000" for red).  Number of provided colors must match the number of
#'               samples in the CAGEset object.  Provided colors are assigned to samples
#'               according to their ordering in the CAGEr object, \emph{i.e} in the order
#'               they are returned by \code{\link{sampleLabels}} function.
#' 
#' @return Assigns one color to each sample in the CAGEr object by setting them as a name
#' attribute of the \code{sampleLabels} slot.
#' 
#' @author Vanja Haberle
#' 
#' @seealso \code{\link{sampleLabels}}
#' @family CAGEr setter methods
#' 
#' @importFrom grDevices rgb
#' 
#' @examples
#' 
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' sampleLabels(exampleCAGEset)
#' setColors(exampleCAGEset, colors = c("darkred", "navy", "forestgreen"))
#' sampleLabels(exampleCAGEset)
#' 
#' @export

setGeneric("setColors", function(object, colors = NULL) standardGeneric("setColors"))

#' @rdname setColors

setMethod("setColors", "CAGEr", function (object, colors){

  objName <- deparse(substitute(object))
  sample.labels <- sampleLabels(object)
  
  if(length(colors) == 1 & is.numeric(colors)){
  	names(sample.labels) <- rainbow(n = length(sample.labels))
  }else if(length(colors) != length(sample.labels)){
  	stop(paste("Number of provided colors must match the number of samples in the CAGEr object, i.e. must be ", length(sample.labels), "!", sep = ""))
  }else if(all(colors %in% colors())){
  	rgb.col <- col2rgb(colors)
  	names(sample.labels) <- apply(rgb.col, 2, function(x) {rgb(red = x[1], green = x[2], blue = x[3], alpha = 255, maxColorValue = 255)})
  }else if((unique(substr(colors, start = 1, stop = 1)) == "#") & all(unique(unlist(strsplit(substr(colors, start = 2, stop = sapply(colors, width)), split = ""))) %in% c(seq(0,9,1), "A", "B", "C", "D", "E", "F"))){
  	names(sample.labels) <- colors
  }else{
  	stop("'colors' argument must be a vector of valid color names in R or a vector of hexadecimal specifications (e.g. #008F0AFF). See colors() for a complete list of valid color names.")
  }
  
  sampleLabels(object) <- sample.labels
  assign(objName, object, envir = parent.frame())
  invisible(1)
})