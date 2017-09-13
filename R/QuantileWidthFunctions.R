#' @name .get.quant.pos
#' 
#' @rdname QuantileWidthFunctions
#' 
#' @title Get quantile positions
#' 
#' @description Function that calculates position of quantiles for CTSS clusters based on
#' distribution of tpm within the cluster
#' 
#' @param cumsums Named list of vectors containing cumulative sum for each cluster
#'        (returned by 'get.cumsum' function).
#' @param clusters GRanges object representing tag clusters or consensus clusters.
#' @param q desired quantiles - single value or a vector of values.
#' @param useMulticore,nrCores See \code{quantilePositions}.
#' 
#' @return Returns the \code{clusters} object with one more metadata column per value
#' in \code{q}, containing Rle integers giving the relative distance of the quantile
#' boundaries to the start position.
#' 
#' @examples 
#' \dontrun{ #because it runs through quantilePositions() anyway.
#' cumsums  <- CTSScumulativesTagClusters(object, 1)
#' clusters <- tagClustersGR(object, 1)
#' .get.quant.pos(cumsums, clusters, c(.1, .9))
#' }

.get.quant.pos <- function(cumsums, clusters, q = NULL, useMulticore = FALSE, nrCores = NULL) {
  getQuantilepos <- function(quant, cumsum) length(Filter(isTRUE, cumsum/max(cumsum) < quant))

	if (.checkMulticore(useMulticore)) {
		cluster.q <- mclapply( cumsums, function(x) sapply(q, getQuantilepos, x)
		                     , mc.cores = .getNrCores(nrCores))
	} else {
		cluster.q <-   lapply(cumsums, function(x) sapply(q, getQuantilepos, x))
	}
	cluster.q <- as.data.frame(do.call(rbind, cluster.q))
	colnames(cluster.q) = paste('q_', q, sep = '')
	mcols(clusters)[,colnames(cluster.q)] <- DataFrame(lapply(cluster.q, Rle))
	clusters
}