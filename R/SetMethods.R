####################################################################
# Functions for setting internal data in CAGEset and CAGEexp objects
#

#' genomeName

setGeneric("genomeName<-", function(object, value) standardGeneric("genomeName<-"))

setMethod("genomeName<-", "CAGEset", function (object, value){
	object@genomeName <- value
	if (validObject(object)) object
})

setMethod("genomeName<-", "CAGEexp", function (object, value){
  metadata(object)$genomeName <- value
	if (validObject(object)) object
})

#' inputFiles

setGeneric("inputFiles<-", function(object, value) standardGeneric("inputFiles<-"))

setMethod("inputFiles<-", "CAGEset", function (object, value){
	object@inputFiles <- value
	if (validObject(object)) object
})

setMethod("inputFiles<-", "CAGEexp", function (object, value){
  object$inputFiles <- value
  if (validObject(object)) object
})

#' inputFilesType

setGeneric("inputFilesType<-", function(object, value) standardGeneric("inputFilesType<-"))

setMethod("inputFilesType<-", "CAGEset", function (object, value){
	object@inputFilesType <- value
	if (validObject(object)) object
})

setMethod("inputFilesType<-", "CAGEexp", function (object, value){
  metadata(object)$inputFilesType <- value
	if (validObject(object)) object
})

#' sampleLabels

setGeneric("sampleLabels<-", function(object, value) standardGeneric("sampleLabels<-"))

setMethod("sampleLabels<-", "CAGEset", function (object, value){
	object@sampleLabels <- value
	if (validObject(object)) object
})

setMethod("sampleLabels<-", "CAGEexp", function (object, value){
  object$sampleLabels <- value
  if (validObject(object)) object
})

#' librarySizes

setGeneric("librarySizes<-", function(object, value) standardGeneric("librarySizes<-"))

setMethod("librarySizes<-", "CAGEset", function (object, value){
	object@librarySizes <- value
	if (validObject(object)) object
})

setMethod("librarySizes<-", "CAGEexp", function (object, value){
  object$librarySizes <- value
  if (validObject(object)) object
})

#' CTSScoordinatesGR

setGeneric("CTSScoordinatesGR<-", function(object, value) standardGeneric("CTSScoordinatesGR<-"))

setMethod("CTSScoordinatesGR<-", "CAGEset", function (object, value){
	stop("Not implemented for the CAGEset class.")
})

setMethod("CTSScoordinatesGR<-", "CAGEexp", function (object, value){
  if (! is(value, "GRanges")) stop("Value must be a GRanges object.")
  rowRanges(object@ExperimentList$tagCountMatrix) <- value
  if (validObject(object)) object
})

#' CTSStagCountSE

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

#' GeneExpSE
#' 
#' Since the SummarizedExperiment can hold normalized and non-normalized values,
#' let's name it "GeneExp" instead of "GeneTagCount" if we would follow the
#' historical CTSS name pattern of CAGEset objects.

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
