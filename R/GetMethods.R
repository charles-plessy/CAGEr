#####################################################
# Functions for retrieving data from CAGEset object
#

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


setGeneric(
name="sampleLabels",
def=function(object){
	standardGeneric("sampleLabels")
})

setMethod("sampleLabels",
signature(object = "CAGEset"),
function (object){
	object@sampleLabels
})

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


setGeneric(
name="tagClusters",
def=function(object, sample){
	standardGeneric("tagClusters")
})

setMethod("tagClusters",
signature(object = "CAGEset"),
function (object, sample){
	if(sample %in% object@sampleLabels){
		object@tagClusters[[sample]]
	}else{
		stop("Provided 'sample' not in the CAGE set! Check sampleLabels()")
	}
})


setGeneric(
name="consensusClusters",
def=function(object){
	standardGeneric("consensusClusters")
})

setMethod("consensusClusters",
signature(object = "CAGEset"),
function (object){
	object@consensusClusters
})


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



###############################################################
# Function for displaying CAGEset object in user friendly way

setMethod("show", 

	signature(object = "CAGEset"), 
	function(object) {
	cat("\nS4 Object of class CAGEset\n\n")
	cat("=======================================\n")
	cat("Input data information\n")
	cat("=======================================\n")
	cat("Reference genome (organism): ", genomeName(object), "\n", sep = "")
	cat("Input file type: ", inputFilesType(object), "\n", sep = "")
	cat("Input file names: ", paste(inputFiles(object), collapse = ", "), "\n", sep = "")
	cat("Sample labels: ", paste(sampleLabels(object), collapse = ", "), "\n", sep = "")
	cat("=======================================\n")
	cat("CTSS information\n")
	cat("=======================================\n")
	if(nrow(CTSScoordinates(object))>0){
	max_out = min(3, nrow(CTSScoordinates(object)))
		cat("CTSS chromosome: ", paste(paste(CTSScoordinates(object)$chr[1:max_out], collapse = ", "), ", ...", sep = ""), "\n", sep = "")
		cat("CTSS position: ", paste(paste(CTSScoordinates(object)$pos[1:max_out], collapse = ", "), ", ...", sep = ""), "\n", sep = "")
		cat("CTSS strand: ", paste(paste(CTSScoordinates(object)$strand[1:max_out], collapse = ", "), ", ...", sep = ""), "\n", sep = "")
	}else{
		cat("CTSS chromosome:\n")	
		cat("CTSS position:\n")	
		cat("CTSS strand:\n")	
	}
	cat("Tag count:\n")
	if(nrow(object@tagCountMatrix) > 0){
		cat(sapply(c(1:length(sampleLabels(object))), function(x) {paste("\t-> ", sampleLabels(object)[x], ": ", paste(paste(object@tagCountMatrix[1:max_out,x], collapse = ", "), ", ...\n", sep = ""), sep = "")}))
	}
	cat("Normalized tpm:\n")
	if(nrow(object@normalizedTpmMatrix)>0){	
		cat(sapply(c(1:length(sampleLabels(object))), function(x) {paste("\t-> ", sampleLabels(object)[x], ": ", paste(paste(formatC(object@normalizedTpmMatrix[1:max_out,x], format = "f", digits = 3), collapse = ", "), ", ...\n", sep = ""), sep = "")}))
	}
	cat("=======================================\n")
	cat("Tag cluster (TC) information\n")
	cat("=======================================\n")
	cat("CTSS clustering method: ", object@clusteringMethod, "\n", sep = "")
	cat("Number of TCs per sample:\n")
	if(length(object@tagClusters) > 0){
		cat(sapply(sampleLabels(object), function(x) {paste("\t-> ", x, ": ", nrow(tagClusters(object, sample = x)), "\n", sep = "")}))
	}
	cat("=======================================\n")
	cat("Consensus cluster information\n")
	cat("=======================================\n")
	if(nrow(consensusClusters(object))>0){
		max_out = min(3, nrow(consensusClusters(object)))
		cat("Number of consensus clusters: ", nrow(consensusClusters(object)), "\n", sep = "")
		cat("Consensus cluster chromosome: ", paste(paste(consensusClusters(object)$chr[1:max_out], collapse = ", "), ", ...", sep = ""), "\n", sep = "")
		cat("Consensus cluster start: ", paste(paste(consensusClusters(object)$start[1:max_out], collapse = ", "), ", ...", sep = ""), "\n", sep = "")
		cat("Consensus cluster end: ", paste(paste(consensusClusters(object)$end[1:max_out], collapse = ", "), ", ...", sep = ""), "\n", sep = "")
		cat("Consensus cluster strand: ", paste(paste(consensusClusters(object)$strand[1:max_out], collapse = ", "), ", ...", sep = ""), "\n", sep = "")
		cat("Normalized tpm:\n")
		cat(sapply(c(1:length(sampleLabels(object))), function(x) {paste("\t-> ", sampleLabels(object)[x], ": ", paste(paste(formatC(object@consensusClustersTpmMatrix[1:max_out,x], format = "f", digits = 3), collapse = ", "), ", ...\n", sep = ""), sep = "")}))
	}else{
		cat("Number of consensus clusters:\n")
		cat("Consensus cluster chromosome:\n")
		cat("Consensus cluster start:\n")
		cat("Consensus cluster end:\n")
		cat("Consensus cluster strand:\n")
		cat("Normalized tpm:\n")		
	}
	cat("=======================================\n")
	cat("Expression profiling\n")
	cat("=======================================\n")
	cat("Expression clustering method: ", object@consensusClustersExpressionClusteringMethod, "\n", sep = "")
	if(length(object@consensusClustersExpressionClasses) > 0){
		cat("Expression clusters for consensus clusters: ", paste(paste(object@consensusClustersExpressionClasses[1:max_out], collapse = ", "), ", ...", sep = ""), "\n", sep = "")
	}else{
		cat("Expression clusters for consensus clusters:\n")
	}
	cat("=======================================\n")
	cat("Promoter shifting\n")
	cat("=======================================\n")
	if(length(object@shiftingGroupX)>0){
		max_out = min(3, nrow(object@consensusClustersShiftingScores))
		cat("GroupX: ", paste(object@shiftingGroupX, collapse = ", "), "\n", sep = "")
		cat("GroupY: ", paste(object@shiftingGroupY, collapse = ", "), "\n", sep = "")
		cat("Shifting scores: ", paste(formatC(object@consensusClustersShiftingScores$shifting.score[1:max_out], format = "f", digits = 3), collapse = ", "), "\n", sep = "")		
	}else{
		cat("GroupX:\n")
		cat("GroupY:\n")
		cat("Shifting scores:\n")
	}	
	cat("\n")
		
})

