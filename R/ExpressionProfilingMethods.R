#' @name getExpressionProfiles
#' 
#' @title CAGE data based expression clustering
#' 
#' @description Clusters CAGE expression across multiple experiments, both at
#' level of individual TSSs or entire clusters of TSSs.
#' 
#' @param object A [`CAGEexp`] object
#' 
#' @param what At which level the expression clustering is done (`CTSS` or
#' `consensusClusters`)
#' 
#' @param tpmThreshold,nrPassThreshold Ignore clusters when their normalized
#' CAGE signal is lower than `tpmThreshold` in at least `nrPassThreshold`
#' experiments.
#' 
#' @param method Method to be used for expression clustering.  `som` uses the
#' self-organizing map (SOM) algorithm of Toronen and coll., FEBS Letters (1999)
#' [`som::som`]] function from _som_ package.  `kmeans` uses the K-means
#' algorithm implemented in the [`stats::kmeans`]] function.
#' 
#' @param xDim,yDim With `method = "kmeans"`, `xDim` specifies number of clusters
#' that will be returned by K-means algorithm and `yDim` is ignored.  With
#' `method = "som"`, `xDim` specifies the the first and `yDim` the second
#' dimension of the self-organizing map, which results in total $xDim x yDim$
#' clusters returned by SOM. 
#' 
#' @details Expression clustering can be done at level of individual CTSSs, in
#' which case the feature vector used as input for clustering algorithm contains
#' log-transformed and scaled (divided by standard deviation) normalized CAGE
#' signal at individual TSS across multiple experiments.  Only TSSs with
#' normalized CAGE signal `>= tpmThreshold` in at least `nrPassThreshold` CAGE
#' experiments are used for expression clustering.  However, CTSSs along the
#' genome can be spatially clustered into tag clusters for each experiment
#' separately using the [`clusterCTSS`] function, and then aggregated across
#' experiments into consensus clusters using [`aggregateTagClusters`] function.
#' Once the consensus clusters have been created, expression clustering at the
#' level of these wider genomic regions (representing entire promoters rather
#' than individual TSSs) can be performed.  In that case the feature vector
#' used as input for clustering algorithm contains normalized CAGE signal
#' within entire consensus cluster across multiple experiments, and threshold
#' values in `tpmThreshold` and `nrPassThreshold` are applied to entire
#' consensus clusters.
#' 
#' @return Returns a modified `CAGEexp` object.  If `what = "CTSS"` the
#' objects's metadata elements `CTSSexpressionClusteringMethod` and
#' `CTSSexpressionClasses` will be set accordingly, and if
#' `what = "consensusClusters"` the elements `consensusClustersExpressionClusteringMethod`
#' and `consensusClustersExpressionClasses` will be set.  Labels of expression
#' classes (clusters) can be retrieved using [`expressionClasses`] function.
#' 
#' @references 
#' Toronen _et al._ (1999) Analysis of gene expression data using
#' self-organizing maps, _FEBS Letters_ *451*:142-146.
#' 
#' @author Vanja Haberle
#' @author Charles Plessy
#' 
#' @family CAGEr expression clustering functions
#' 
#' @examples 
#' getExpressionProfiles( exampleCAGEexp, "CTSS"
#'                      , tpmThreshold = 50, nrPassThreshold = 1
#'                      , method = "som", xDim = 3, yDim = 3)
#'                      
#' getExpressionProfiles( exampleCAGEexp, "CTSS"
#'                      , tpmThreshold = 50, nrPassThreshold = 1
#'                      , method = "kmeans", xDim = 3)
#' 
#' getExpressionProfiles(exampleCAGEexp, "consensusClusters")
#' 
#' @importFrom stats kmeans
#' @export

setGeneric( "getExpressionProfiles"
          , function( object, what = c("CTSS", "consensusClusters")
                    , tpmThreshold = 5, nrPassThreshold = 1
                    , method = c("som", "kmeans")
                    , xDim = 5, yDim = 5)
              standardGeneric("getExpressionProfiles"))

#' @rdname getExpressionProfiles

setMethod( "getExpressionProfiles", "CAGEexp"
         , function (object, what, tpmThreshold, nrPassThreshold, method, xDim, yDim){
  what <- match.arg(what)

	if (length(sampleLabels(object)) < 2)
		stop("Provided CAGEexp object contains only one sample! At least two samples are required for expression profiling!")
  
  tpm.mx <- switch( what
                  , CTSS              = CTSSnormalizedTpmDF(object)
                  , consensusClusters = assay(consensusClustersSE(object), "normalized"))
  
  l   <- getExpressionProfiles( object = DelayedArray(tpm.mx)
                              , tpmThreshold = tpmThreshold, nrPassThreshold = nrPassThreshold
                              , method = method, xDim = xDim, yDim = yDim)
  cl  <- l[[1]]
  idx <- l[[2]]

  if (what == "CTSS") {
    mcols(CTSScoordinatesGR(object))$exprClass <- Rle(NA) # initialise
    mcols(CTSScoordinatesGR(object))[idx, "exprClass"] <- cl
    metadata(object)$CTSSexpressionClusteringMethod <- method	
 } else if (what == "consensusClusters") {
    mcols(consensusClustersGR(object))$exprClass <- Rle(NA) # initialise
    mcols(consensusClustersGR(object))[idx, "exprClass"] <- cl
    metadata(object)$consensusClustersExpressionClusteringMethod <- method
  }
  object
})

#' @rdname getExpressionProfiles
#' @importFrom som som
#' @importFrom stats kmeans
#' @import DelayedMatrixStats

setMethod("getExpressionProfiles", "DelayedArray",
function (object, what, tpmThreshold, nrPassThreshold, method, xDim, yDim) {
  method <- match.arg(method)

  idx <- .filterCtss( object, threshold = tpmThreshold
                    , nrPassThreshold = nrPassThreshold, thresholdIsTpm = TRUE)
	
	m <- t(scale(t(log(object + 1)), center=FALSE))
	m <- m[idx,]
	
	if(method == "som") {
		m.som <- som(m, xDim, yDim, neigh = "gaussian", topol = "hex")
		cl <- paste(m.som$visual$x, m.som$visual$y, sep = "_")
	}else if(method == "kmeans"){
		m.kmeans <- kmeans(m, centers = xDim)
		cl <- as.character(m.kmeans$cluster)
	}
	return(list(cl, idx))
})