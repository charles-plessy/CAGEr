setGeneric(
name="clusterCTSS",
def=function(object, threshold = 1, nrPassThreshold = 1, thresholdIsTpm = TRUE, method = "distclu", maxDist = 20, removeSingletons = FALSE, keepSingletonsAbove = Inf, minStability = 1, maxLength = 500, reduceToNonoverlapping = TRUE, customClusters = NULL, useMulticore = FALSE, nrCores = NULL){
	standardGeneric("clusterCTSS")
}
)

setMethod("clusterCTSS",
signature(object = "CAGEset"),
function (object, threshold = 1, nrPassThreshold = 1, thresholdIsTpm = TRUE, method = "distclu", maxDist = 20, removeSingletons = FALSE, keepSingletonsAbove = Inf, minStability = 1, maxLength = 500, reduceToNonoverlapping = TRUE, customClusters = NULL, useMulticore = FALSE, nrCores = NULL){
	
	pt <- .Platform$OS.type
	if(useMulticore == TRUE){
		if(pt == "unix"){
			if("parallel" %in% rownames(installed.packages()) == FALSE){
				stop("Cannot use multicore because package 'parallel' is not installed!")
			}
		}else{
			useMulticore = FALSE
			warning("Multicore is not supported on non-Unix platforms! Setting useMulticore=FALSE")
		}
	}
		
	objName <- deparse(substitute(object))
	sample.labels <- sampleLabels(object)
	
	message("\nFiltering CTSSs below threshold...")
	if(thresholdIsTpm){
		data <- object@normalizedTpmMatrix
	}else{
		data <- object@tagCountMatrix
	}
	
	if(threshold > 0){
		nr.pass.threshold <- apply(data, 1, function(x) {sum(x >= threshold)})
		idx <- nr.pass.threshold >= min(nrPassThreshold, length(sample.labels))
		gc()
		data <- cbind(CTSScoordinates(object)[idx,], object@normalizedTpmMatrix[idx,,drop=F])
		gc()
	}else{
		data <- cbind(CTSScoordinates(object), object@normalizedTpmMatrix)
		idx <- rep(TRUE, nrow(data))

	}
	
	message("Clustering...")
	if(method == "distclu"){
		ctss.cluster.list <- .distclu(data = data, sample.labels = sample.labels, max.dist = maxDist, removeSingletons = removeSingletons, keepSingletonsAbove = keepSingletonsAbove, useMulticore = useMulticore, nrCores = nrCores)
	}else if (method == "paraclu"){
		ctss.cluster.list <- .paraclu(data = data, sample.labels = sample.labels, minStability = minStability, maxLength = maxLength, removeSingletons = removeSingletons, keepSingletonsAbove = keepSingletonsAbove, reduceToNonoverlapping = reduceToNonoverlapping, useMulticore = useMulticore, nrCores = nrCores)
	}else if(method == "custom"){
		if(length(customClusters)==0){
			stop("'customClusters' must be given when method = \"custom\"")
		}
		ctss.cluster.list <- .predefined.clusters(data = data, sample.labels = sample.labels, custom.clusters = customClusters, useMulticore = useMulticore, nrCores = nrCores)
	}else{
		stop("'method' parameter must be one of the (\"distclu\", \"paraclu\", \"custom\")")
	}
	
	object@filteredCTSSidx <- idx
	object@clusteringMethod <- method
	object@tagClusters <- ctss.cluster.list
	assign(objName, object, envir = parent.frame())
	invisible(1)
	
}
)

