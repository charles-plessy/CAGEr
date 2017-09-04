#' @include CAGEr.R

#' @name quantilePositions
#' 
#' @title Determining positions of CAGE signal quantiles within genomic region
#' 
#' @description Calculates positions of quantiles of CAGE signal along tag clusters or consensus
#' clusters in each CAGE dataset within CAGEset object.  The function calculates
#' positions of both \dQuote{lower} and \dQuote{upper} quantiles as described in Details.
#' 
#' @param object 	A \code{\link{CAGEset}} object
#' 
#' @param clusters 	Which clusters should be used.  Can be either \code{tagClusters}
#' to calculate positions of quantiles in tag clusters (different set of genomic
#' coordinates for every CAGE experiment) or \code{consensusClusters} to calculate
#' positions of quantiles in consensus clusters (same set of genomic coordinates for
#' every CAGE experiment).
#' 
#' @param qLow Which "lower" quantiles should be calculated.  It has to be a numeric
#' vector of values in range \code{[0,1]}.  See Details.
#' 
#' @param qUp Which "upper" quantiles should be calculated.  It has to be a numeric
#' vector of values in range \code{[0,1]}.  See Details.
#' 
#' @param useMulticore Logical, should multicore be used.  \code{useMulticore = TRUE}
#' is supported only on Unix-like platforms.
#' 
#' @param nrCores Number of cores to use when \code{useMulticore = TRUE}.  Default
#' value \code{NULL} uses all detected cores.
#' 
#' @details Position of the "lower" quantile \code{qLow} is defined as a point that divides the
#' genomic region into two parts, so that the 5' part contains \code{< qLow * 100\%} of
#' the CAGE signal of that region. Accordingly, position of the "upper" quantile
#' \code{qUp} is defined as a point that divides the genomic region into two parts so
#' that the 5' part contains \code{>= qUp * 100\%} of the CAGE signal of that region.
#' Positions of one "lower" and one "upper" quantile (when \code{qLow <= qUp}) define a
#' central part of the genomic region that contains \code{>= (qUp - qLow) * 100\%} of
#' the CAGE signal of that region. Width of that central part is refered to as
#' "interquantile width", which is a more robust measure of the promoter width than the
#' total span of the region.
#' 
#' @return When \code{clusters = "tagClusters"}, the slots \code{tagClustersQuantileLow}
#' and \code{tagClustersQuantileUp} of the provided \code{\link{CAGEset}} object will
#' be occupied with the positions of specified quantiles in all tag clusters for all
#' CAGE datasets. When \code{clusters = "consensusClusters"} the slots
#' \code{consensusClustersQuantileLow} and \code{consensusClustersQuantileUp} will be
#' occupied by the corresponding information for consensus clusters.
#' 
#' @author Vanja Haberle
#' 
#' @family CAGEr object modifiers
#' @family CAGEr clusters functions
#' 
#' @examples 
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' head(cbind(
#'   CAGEr:::tagClustersQuantileLow(exampleCAGEset)[[1]],
#'   CAGEr:::tagClustersQuantileUp(exampleCAGEset)[[1]]
#' ))
#' quantilePositions( object = exampleCAGEset, clusters = "tagClusters"
#'                  , qLow = c(0.1,0.2), qUp = c(0.8,0.9))
#' head(cbind(
#'   CAGEr:::tagClustersQuantileLow(exampleCAGEset)[[1]],
#'   CAGEr:::tagClustersQuantileUp(exampleCAGEset)[[1]]
#' ))
#' 
#' quantilePositions( object = exampleCAGEset, clusters = "consensusClusters"
#'                  , qLow = c(0.1,0.2), qUp = c(0.8,0.9))
#'                  
#' head(cbind(
#'   CAGEr:::consensusClustersQuantileLow(exampleCAGEset)[[1]],
#'   CAGEr:::consensusClustersQuantileUp(exampleCAGEset)[[1]]
#' ))
#' 
#' ce <- readRDS(system.file(package = "CAGEr", "extdata/CAGEexp.rds"))
#' normalizeTagCount(ce)
#' clusterCTSS( object = ce, threshold = 50, thresholdIsTpm = TRUE
#'            , nrPassThreshold = 1, method = "distclu", maxDist = 20
#'            , removeSingletons = TRUE, keepSingletonsAbove = 100)
#' cumulativeCTSSdistribution(ce, clusters = "tagClusters")
#' quantilePositions( ce, clusters = "tagClusters"
#'                  , qLow = c(0.1,0.2), qUp = c(0.8,0.9))
#' 
#' @export

setGeneric( "quantilePositions"
          , function(object, clusters, qLow = 0.1, qUp = 0.9, useMulticore = FALSE, nrCores = NULL)
	standardGeneric("quantilePositions"))

#' @rdname quantilePositions

setMethod("quantilePositions",
signature(object = "CAGEr"),
function (object, clusters, qLow, qUp, useMulticore, nrCores){

  useMulticore <- .checkMulticore(useMulticore)
	
	objName <- deparse(substitute(object))

	message("\nGetting positions of quantiles within clusters...")
	
	ctss.clusters.q.low.list <- list()
	ctss.clusters.q.up.list <- list()
	
	if(clusters == "tagClusters"){	
		
		samples.cumsum.list <- CTSScumulativesTagClusters(object)
		
		for(s in sampleLabels(object)) {
			message("\t-> ", s)
			
			clusters.cumsum.list <- samples.cumsum.list[[s]]
			ctss.clusters <- tagClusters(object, sample = s)
			ctss.clusters.q.low <- .get.quant.pos(cluster.cumsums = clusters.cumsum.list, coors = ctss.clusters, q = qLow, useMulticore = useMulticore, nrCores = nrCores)
			ctss.clusters.q.up <- .get.quant.pos(cluster.cumsums = clusters.cumsum.list, coors = ctss.clusters, q = qUp, useMulticore = useMulticore, nrCores = nrCores)

			ctss.clusters.q.low.list[[s]] <- ctss.clusters.q.low[,c(which(colnames(ctss.clusters.q.low) == "cluster"), grep("q_", colnames(ctss.clusters.q.low), fixed = T))]
			ctss.clusters.q.up.list[[s]] <- ctss.clusters.q.up[,c(which(colnames(ctss.clusters.q.up) == "cluster"), grep("q_", colnames(ctss.clusters.q.up), fixed = T))]
		
		}
		
		tagClustersQuantileLow(object) <- ctss.clusters.q.low.list
		tagClustersQuantileUp(object)  <- ctss.clusters.q.up.list
		
	}else if (clusters == "consensusClusters"){
		samples.cumsum.list <- CTSScumulativesCC(object)
		ctss.clusters.orig <- merge(consensusClusters(object), consensusClustersTpm(object), by.x = 1, by.y = 0)
		
		for(s in sampleLabels(object)) {
			message("\t-> ", s)
			
			clusters.cumsum.list <- samples.cumsum.list[[s]]
			ctss.clusters <- ctss.clusters.orig[ctss.clusters.orig[,s]>0,]
			colnames(ctss.clusters)[which(colnames(ctss.clusters) == "consensus.cluster")] = "cluster"
			ctss.clusters.q.low <- .get.quant.pos(cluster.cumsums = clusters.cumsum.list, coors = ctss.clusters, q = qLow, useMulticore = useMulticore, nrCores = nrCores)
			ctss.clusters.q.up <- .get.quant.pos(cluster.cumsums = clusters.cumsum.list, coors = ctss.clusters, q = qUp, useMulticore = useMulticore, nrCores = nrCores)
			
			ctss.clusters.q.low.list[[s]] <- ctss.clusters.q.low[,c(which(colnames(ctss.clusters.q.low) == "cluster"), grep("q_", colnames(ctss.clusters.q.low), fixed = T))]
			ctss.clusters.q.up.list[[s]] <- ctss.clusters.q.up[,c(which(colnames(ctss.clusters.q.up) == "cluster"), grep("q_", colnames(ctss.clusters.q.up), fixed = T))]
		}
		consensusClustersQuantileLow(object) <- ctss.clusters.q.low.list
		consensusClustersQuantileUp(object)  <- ctss.clusters.q.up.list
	}else{
		stop("'clusters' parameter must be one of the (\"tagClusters\", \"consensusClusters\")")
	}
	assign(objName, object, envir = parent.frame())
	invisible(1)
})
