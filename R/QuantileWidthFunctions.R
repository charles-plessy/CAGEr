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

.get.quant.pos <- function(cumsums, clusters, q) {
  getQuantilepos <- Vectorize(vectorize.args = "cumsum", function(q, cumsum) {
    cumsum <- decode(cumsum)
    max <- tail(cumsum,1)   # Max is last element since we x is a cumulative sums.
    treshold <- max * q
    which.max(cumsum >= treshold)
  })
	cluster.q <- lapply(q, getQuantilepos, cumsums)
	names(cluster.q) = paste('q_', q, sep = '')
	mcols(clusters)[, names(cluster.q)] <- DataFrame(lapply(cluster.q, Rle))
	clusters
}