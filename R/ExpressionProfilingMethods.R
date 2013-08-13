###################################################################
# Functions for clustering by expression (expression profiling)
#

setGeneric(
name="getExpressionProfiles",
def=function(object, what, tpmThreshold = 5, nrPassThreshold = 1, method = "som", xDim = 5, yDim = 5){
	standardGeneric("getExpressionProfiles")
}
)

setMethod("getExpressionProfiles",
signature(object = "CAGEset"),
function (object, what, tpmThreshold, nrPassThreshold, method, xDim, yDim){

	objName <- deparse(substitute(object))
	sample.labels = sampleLabels(object)
	
	if(length(sample.labels) < 2){
		stop("Provided CAGEset object contains only one sample! At least two samples are required for expression profiling!")
	}

	if(what == "CTSS") {
		
		tpm.mx <- object@normalizedTpmMatrix
		l = .clusterExpression(tpm.mx = tpm.mx, sample.labels = sample.labels, tpmThreshold = tpmThreshold, nrPassThreshold = nrPassThreshold, method = method, xDim = xDim, yDim = yDim)
		cl = l[[1]]
		idx = l[[2]]
		names(cl) <- which(idx == TRUE)
		object@CTSSexpressionClasses <- cl
		object@CTSSexpressionClusteringMethod <- method		
	
	}else if(what == "consensusClusters") {
	
		tpm.mx <- object@consensusClustersTpmMatrix
		l = .clusterExpression(tpm.mx = tpm.mx, sample.labels = sample.labels, tpmThreshold = tpmThreshold, nrPassThreshold = nrPassThreshold, method = method, xDim = xDim, yDim = yDim)
		cl = l[[1]]
		idx = l[[2]]
		names(cl) <- which(idx == TRUE)
		object@consensusClustersExpressionClasses <- cl
		object@consensusClustersExpressionClusteringMethod <- method
		
	}else{
		stop("'what' parameter must be one of the (\"CTSS\", \"consensusClusters\")")
	}

	assign(objName, object, envir = parent.frame())
	invisible(1)	
		
}
)


.clusterExpression <- function(tpm.mx, sample.labels, tpmThreshold = 5, nrPassThreshold = 1, method, xDim = 5, yDim = 5) {
	
	nr.pass.threshold <- apply(tpm.mx, 1, function(x) {sum(x >= tpmThreshold)})
	idx <- nr.pass.threshold >= min(nrPassThreshold, length(sample.labels))
	
	m <- t(scale(t(log(tpm.mx+1)), center=F)) 
	m <- m[idx,]
	
	if(method == "som") {
		m.som <- som(m, xDim, yDim, neigh = "gaussian", topol = "hex")
		cl <- paste(m.som$visual$x, m.som$visual$y, sep = "_")
	}else if(method == "kmeans"){
		m.kmeans <- kmeans(m, centers = xDim)
		cl <- as.character(m.kmeans$cluster)
	}else{
		stop("'method' parameter must be one of the (\"som\", \"kmeans\")")
	}	
	return(list(cl, idx))
	
}



setGeneric(
name="extractExpressionClass",
def=function(object, what, which="all"){
	standardGeneric("extractExpressionClass")
})

setMethod("extractExpressionClass",
signature(object = "CAGEset"),
function (object, what, which="all"){
	
	objName <- deparse(substitute(object))

	if(what == "CTSS"){
		classes <- object@CTSSexpressionClasses
		if(length(classes)>0){
			r <- cbind(CTSScoordinates(object), object@normalizedTpmMatrix)
			r$expression_class <- NA
			r$expression_class[as.integer(names(classes))] <- classes
			if(which == "all") {
				r <- subset(r, !(is.na(expression_class)))
				return(r)
			}else if(as.character(which) %in% unique(classes)){
				r <- subset(r, expression_class == as.character(which))
				return(r)
			}else{
				stop("'which' parameter must be either 'all' or one of the expression classes! See expressionClasses(", objName, ", what='CTSS') for CTSS expression classes in your dataset!")
			}
		}else{
			stop("No CTSS expression clustering has been done yet!")
		}
		
	}else if(what == "consensusClusters") {
		classes <- object@consensusClustersExpressionClasses
		if(length(classes)>0){
			r <- cbind(consensusClusters(object),  object@consensusClustersTpmMatrix)
			r$expression_class <- NA
			r$expression_class[as.integer(names(classes))] <- classes
			if(which == "all") {
				r <- subset(r, !(is.na(expression_class)))
				return(r)
			}else if(as.character(which) %in% unique(classes)){
				r <- subset(r, expression_class == as.character(which))
				return(r)
			}else{
				stop("'which' parameter must be either 'all' or one of the expression classes! See expressionClasses(", objName, ", what='consensusClusters') for consensusCluster expression classes in your dataset!")
			}
		}else{
			stop("No consensusClusters expression clustering has been done yet!")
		}

	}else{
		stop("'what' parameter must be one of the (\"CTSS\", \"consensusClusters\")")
	}

})




















