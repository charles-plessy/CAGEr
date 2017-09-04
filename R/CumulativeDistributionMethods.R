#' @include GetMethods.R

#' @name cumulativeCTSSdistribution
#' 
#' @title Calculating cumulative sum of CAGE signal along genomic region
#' 
#' @description Calculates cumulative sum of CAGE signal along each tag cluster or
#' consensus cluster in every CAGE dataset within CAGEset object.
#' 
#' @param object A \code{\link{CAGEset}} object
#' @param clusters Which clusters should be used.  Can be either \code{clusters = "tagClusters"}
#'        to calculate cumulative sum along tag clusters (different set of genomic coordinates
#'        for every CAGE experiment) or \code{clusters = "consensusClusters"} to calculate
#'        cumulative sum along consensus clusters (same set of genomic coordinates for every
#'        CAGE experiment).
#' @param useMulticore Logical, should multicore be used.  \code{useMulticore = TRUE} is
#'        supported only on Unix-like platforms.
#' @param nrCores Number of cores to use when \code{useMulticore = TRUE}.  Default value
#'        \code{NULL} uses all detected cores.
#' 
#' @return The slot \code{CTSScumulativesTagClusters} (when \code{clusters = "tagClusters"}) or
#' \code{CTSScumulativesConsensusClusters} (when \code{clusters = "consensusClusters"}) of the
#' provided \code{\link{CAGEset}} object will be occupied by the list containing cumulative sum
#' of the CAGE signal along genomic regions per CAGE experiment.
#' 
#' @author Vanja Haberle
#' 
#' @family CAGEr object modifiers
#' @family CAGEr clusters functions
#' 
#' @examples
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' CTSScumulativesTagClusters(exampleCAGEset)[[1]][1:6]
#' cumulativeCTSSdistribution(object = exampleCAGEset, clusters = "tagClusters")
#' CTSScumulativesTagClusters(exampleCAGEset)[[1]][1:6]
#' cumulativeCTSSdistribution(exampleCAGEset, clusters = "consensusClusters")
#' 
#' ce <- readRDS(system.file(package = "CAGEr", "extdata/CAGEexp.rds"))
#' normalizeTagCount(ce)
#' clusterCTSS( object = ce, threshold = 50, thresholdIsTpm = TRUE
#'            , nrPassThreshold = 1, method = "distclu", maxDist = 20
#'            , removeSingletons = TRUE, keepSingletonsAbove = 100)
#' cumulativeCTSSdistribution(ce, clusters = "tagClusters")
#' cumulativeCTSSdistribution(ce, clusters = "consensusClusters")
#'            
#' @export

setGeneric( "cumulativeCTSSdistribution"
          , function ( object, clusters
                     , useMulticore = FALSE, nrCores = NULL) {
              standardGeneric("cumulativeCTSSdistribution")})

setMethod("cumulativeCTSSdistribution", "CAGEr",
function (object, clusters, useMulticore = FALSE, nrCores = NULL){
  objName <- deparse(substitute(object))
	useMulticore <- .checkMulticore(useMulticore)
	message("\nCalculating cumulative sum of CAGE signal along clusters...")
	samples.cumsum.list <- list()
	if(clusters == "tagClusters"){	
		for(s in sampleLabels(object)) {
			message("\t-> ", s)
			samples.cumsum.list[[s]] <-
			  .getCumsum( ctss      = .CTSS(CTSSnormalizedTpmGR(object, s))
			            , clusters  = getTagClusterGR(object, s)
			            , useMulticore = useMulticore, nrCores = nrCores)
		}
		CTSScumulativesTagClusters(object) <-samples.cumsum.list
	}else if (clusters == "consensusClusters"){
		for(s in sampleLabels(object)) {
			message("\t-> ", s)
		  	samples.cumsum.list[[s]] <-
			  .getCumsum( ctss      = .CTSS(CTSSnormalizedTpmGR(object, s))
			            , clusters  = consensusClustersGR(object)
			            , useMulticore = useMulticore, nrCores = nrCores)
		}
		CTSScumulativesCC(object) <-samples.cumsum.list
		
	}else{
		stop("'clusters' parameter must be one of the (\"tagClusters\", \"consensusClusters\")")
	}
	assign(objName, object, envir = parent.frame())
	invisible(1)
})
