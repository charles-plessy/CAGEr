####################################################################
# Functions for setting internal data in CAGEset and CAGEexp objects
#

#' @name genomeName
#' @aliases genomeName<-
#' @family CAGEr setter methods
#' @docType methods
#' @author Charles Plessy
#' @export

setGeneric("genomeName<-", function(object, value) standardGeneric("genomeName<-"))

setMethod("genomeName<-", "CAGEset", function (object, value){
	object@genomeName <- value
	if (validObject(object)) object
})

setMethod("genomeName<-", "CAGEexp", function (object, value){
  metadata(object)$genomeName <- value
	if (validObject(object)) object
})

#' @name inputFiles
#' @aliases inputFiles<-
#' @family CAGEr setter methods
#' @docType methods
#' @author Charles Plessy
#' @export

setGeneric("inputFiles<-", function(object, value) standardGeneric("inputFiles<-"))

setMethod("inputFiles<-", "CAGEset", function (object, value){
	object@inputFiles <- value
	if (validObject(object)) object
})

setMethod("inputFiles<-", "CAGEexp", function (object, value){
  object$inputFiles <- value
  if (validObject(object)) object
})

#' @name inputFilesType
#' @aliases inputFilesType<-
#' @family CAGEr setter methods
#' @docType methods
#' @author Charles Plessy
#' @export

setGeneric("inputFilesType<-", function(object, value) standardGeneric("inputFilesType<-"))

setMethod("inputFilesType<-", "CAGEset", function (object, value){
	object@inputFilesType <- value
	if (validObject(object)) object
})

setMethod("inputFilesType<-", "CAGEexp", function (object, value){
  object$inputFilesType <- value
	if (validObject(object)) object
})

#' @name sampleLabels
#' @aliases sampleLabels<-
#' @family CAGEr setter methods
#' @docType methods
#' @author Charles Plessy
#' @export

setGeneric("sampleLabels<-", function(object, value) standardGeneric("sampleLabels<-"))

setMethod("sampleLabels<-", "CAGEset", function (object, value){
	object@sampleLabels <- value
	if (validObject(object)) object
})

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

#' @name CTSScoordinatesGR
#' @noRd
#' @export

setGeneric("CTSScoordinatesGR<-", function(object, value) standardGeneric("CTSScoordinatesGR<-"))

setMethod("CTSScoordinatesGR<-", "CAGEset", function (object, value){
	stop("Not implemented for the CAGEset class.")
})

setMethod("CTSScoordinatesGR<-", "CAGEexp", function (object, value){
  if (! is(value, "GRanges")) stop("Value must be a GRanges object.")
  rowRanges(object@ExperimentList$tagCountMatrix) <- value
  if (validObject(object)) object
})

#' @name CTSStagCountSE
#' @noRd
#' @export

setGeneric("CTSStagCountSE<-", function(object, value) standardGeneric("CTSStagCountSE<-"))

setMethod("CTSStagCountSE<-", "CAGEset", function (object, value){
	stop("Not implemented for the CAGEset class.")
})

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


#' @name `CTSScumulativesTagClusters<-`
#' 
#' @rdname CTSScumulativesTagClusters
#' 
#' @param object A \code{\link{CAGEset}} or \code{\link{CAGEset}} object.
#' @param value CTSScumulativesTagClusters data

setGeneric( "CTSScumulativesTagClusters<-"
          , function(object, value) standardGeneric("CTSScumulativesTagClusters<-"))

setMethod("CTSScumulativesTagClusters<-", "CAGEset", function (object, value) {
	object@CTSScumulativesTagClusters <- value
	if (validObject(object)) object
})

setMethod("CTSScumulativesTagClusters<-", "CAGEexp", function (object, value) {
  metadata(object)$CTSScumulativesTagClusters <- value
  if (validObject(object)) object
})

#' @name `tagClustersQuantileLow<-`
#' @rdname tagClustersQuantile
#' 

setGeneric("tagClustersQuantileLow<-", function(object, value) standardGeneric("tagClustersQuantileLow<-"))

setMethod("tagClustersQuantileLow<-", "CAGEset", function (object, value){
	object@tagClustersQuantileLow <- value
	if (validObject(object)) object
})

setMethod("tagClustersQuantileLow<-", "CAGEexp", function (object, value){
  metadata(object)$tagClustersQuantileLow <- value
  if (validObject(object)) object
})

#' @name tagClustersQuantileUp
#' @rdname tagClustersQuantile
#' 

setGeneric("tagClustersQuantileUp<-", function(object, value) standardGeneric("tagClustersQuantileUp<-"))

setMethod("tagClustersQuantileUp<-", "CAGEset", function (object, value){
	object@tagClustersQuantileUp <- value
  if (validObject(object)) object
})

setMethod("tagClustersQuantileUp<-", "CAGEexp", function (object, value){
  metadata(object)$tagClustersQuantileUp <- value
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
#' @examples
#' 
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' sampleLabels(exampleCAGEset)
#' setColors(exampleCAGEset, colors = c("darkred", "navy", "forestgreen"))
#' sampleLabels(exampleCAGEset)
#' 
#' @export

setGeneric(
name="setColors",
def=function(object, colors = NULL){
	standardGeneric("setColors")
}
)

setMethod("setColors",
signature(object = "CAGEr"),
function (object, colors = NULL){

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