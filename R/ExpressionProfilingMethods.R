#' @name getExpressionProfiles
#' 
#' @title CAGE data based expression clustering
#' 
#' @description Clusters CAGE expression across multiple experiments, both at level of
#' individual TSSs or entire clusters of TSSs.
#' 
#' @param object A \code{\link{CAGEexp}} object
#' 
#' @param what At which level should the expression clustering be done.  Can be either
#' \code{"CTSS"} to perform clustering of individual CTSSs or \code{"consensusClusters"}
#' to perform clustering of consensus clusters.  See Details.
#' 
#' @param tpmThreshold,nrPassThreshold Only CTSSs or consensus clusters
#' (depending on \code{what} parameter) with normalized CAGE signal \code{>= tpmThreshold}
#' in \code{>= nrPassThreshold} experiments will be included in expression clustering.
#' 
#' @param method Method to be used for expression clustering.  Can be either \code{"som"} to
#' use the self-organizing map (SOM) algorithm (Toronen \emph{et al}., FEBS Letters 1999)
#' implemented in the the \code{som} function from \code{som} package, or \code{"kmeans"}
#' to use the K-means algorithm implemented in the \code{kmeans} function from \code{stats}
#' package.
#' 
#' @param xDim,yDim When \code{method = "kmeans"}, \code{xDim} specifies number of clusters
#' that will be returned by K-means algorithm and \code{yDim} is ignored.  When \code{method =
#' "som"}, \code{xDim} specifies the the first and \code{yDim} the second dimension of the self
#' -organizing map, which results in total \code{xDim * yDim} clusters returned by SOM.
#' 
#' @details Expression clustering can be done at level of individual CTSSs, in which case the
#' feature vector used as input for clustering algorithm contains log-transformed and scaled
#' (divided by standard deviation) normalized CAGE signal at individual TSS across multiple
#' experiments.  Only TSSs with normalized CAGE signal \code{>= tpmThreshold} in at least
#' \code{nrPassThreshold} CAGE experiments are used for expression clustering.  However, CTSSs
#' along the genome can be spatially clustered into tag clusters for each experiment separately
#' using the \code{\link{clusterCTSS}} function, and then aggregated across experiments into
#' consensus clusters using \code{\link{aggregateTagClusters}} function.  Once the consensus
#' clusters have been created, expression clustering at the level of these wider genomic regions
#' (representing entire promoters rather than individual TSSs) can be performed.  In that case
#' the feature vector used as input for clustering algorithm contains normalized CAGE signal
#' within entire consensus cluster across multiple experiments, and threshold values in
#' \code{tpmThreshold} and \code{nrPassThreshold} are applied to entire consensus clusters.
#' 
#' @return If \code{what = "CTSS"} the slots \code{CTSSexpressionClusteringMethod} and
#' \code{CTSSexpressionClasses} will be occupied, and if \code{what = "consensusClusters"} the
#' slots \code{consensusClustersExpressionClusteringMethod} and
#' \code{consensusClustersExpressionClasses} of the provided \code{\link{CAGEexp}} object will
#' be occupied with the results of expression clustering.  Labels of expression classes (clusters)
#' can be retrieved using \code{\link{expressionClasses}} function, and elements belonging to a
#' specific expression class can be selected using \code{\link{extractExpressionClass}} function.
#' 
#' @references 
#' Toronen \emph{et al}. (1999) Analysis of gene expression data using self-organizing maps,
#' \emph{FEBS Letters} \bold{451}:142-146.
#' 
#' @author Vanja Haberle
#' 
#' @seealso \code{\link{plotExpressionProfiles}}, \code{\link{expressionClasses}},
#' \code{\link{extractExpressionClass}}.
#' 
#' @examples 
#' getExpressionProfiles(exampleCAGEexp, what = "CTSS", tpmThreshold = 50, nrPassThreshold = 1
#'                      , method = "som", xDim = 3, yDim = 3)
#' 
#' @importFrom stats kmeans
#' @export

setGeneric( "getExpressionProfiles"
          , function( object, what
                    , tpmThreshold = 5, nrPassThreshold = 1
                    , method = "som"
                    , xDim = 5, yDim = 5)
              standardGeneric("getExpressionProfiles"))

#' @rdname getExpressionProfiles

setMethod( "getExpressionProfiles", "CAGEexp"
         , function (object, what, tpmThreshold, nrPassThreshold, method, xDim, yDim){

	objName <- deparse(substitute(object))
	sample.labels = sampleLabels(object)
	
	if(length(sample.labels) < 2){
		stop("Provided CAGEexp object contains only one sample! At least two samples are required for expression profiling!")
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
})

#' @name .clusterExpression
#' @noRd
#' @importFrom som som

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

#' @name extractExpressionClass
#' 
#' @title Extract elements of the specified expression class
#' 
#' @description Extracts CTSSs or consensus clusters belonging to a specified expression class.
#' 
#' @param object A \code{\link{CAGEexp}} object.
#' 
#' @param what Which level of expression clustering should be used. Can be either
#'        \code{"CTSS"} to extract expression class of individual CTSSs or
#'        \code{"consensusClusters"} to extract expression class of consensus clusters.
#' 
#' @param which Which expression class should be extracted. It has to be one of the valid
#'        expression class labels (as returned by \code{\link{expressionClasses}} function),
#'        or \code{"all"} to extract members of all expression classes.
#'        
#' @return Returns a \code{data.frame} of CTSSs (when \code{what = "CTSS"}) or consensus clusters
#' (when \code{what = "consensusClusters"}) belonging to a specified expression class, with
#' genomic coordinates and additional information.  Last column contains the label of the
#' corresponding expression class.
#' 
#' @author Vanja Haberle
#' 
#' @seealso \code{\link{getExpressionProfiles}}, \code{\link{plotExpressionProfiles}},
#'          \code{\link{expressionClasses}}.
#' 
#' @examples
#' CTSSexprClasses <- extractExpressionClass(exampleCAGEexp, what = "CTSS", which = "all")
#' head(CTSSexprClasses)
#' 
#' @export

setGeneric( "extractExpressionClass", function(object, what, which="all")
  standardGeneric("extractExpressionClass"))

#' @rdname extractExpressionClass

setMethod( "extractExpressionClass", "CAGEexp", function (object, what, which="all")
  stop("Not supported for CAGEexp objects."))

#' @rdname extractExpressionClass

setMethod( "extractExpressionClass", "CAGEr", function (object, what, which="all"){
	objName <- deparse(substitute(object))

	if(what == "CTSS"){
		classes <- object@CTSSexpressionClasses
		if(length(classes)>0){
			r <- cbind(CTSScoordinates(object), object@normalizedTpmMatrix)
			r$expression_class <- NA
			r$expression_class[as.integer(names(classes))] <- classes
			if(which == "all") {
				r <- subset(r, !(is.na(r$expression_class)))
				return(r)
			}else if(as.character(which) %in% unique(classes)){
				r <- subset(r, r$expression_class == as.character(which))
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
				r <- subset(r, !(is.na(r$expression_class)))
				return(r)
			}else if(as.character(which) %in% unique(classes)){
				r <- subset(r, r$expression_class == as.character(which))
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
