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
#' head(cbind(
#'   CAGEr:::tagClustersQuantileLow(exampleCAGEset, 1),
#'   CAGEr:::tagClustersQuantileUp (exampleCAGEset, 1)
#' ))
#' quantilePositions( object = exampleCAGEset, clusters = "tagClusters"
#'                  , qLow = c(0.1, 0.2), qUp = c(0.8, 0.9))
#' head(cbind(
#'   CAGEr:::tagClustersQuantileLow(exampleCAGEset, 1),
#'   CAGEr:::tagClustersQuantileUp (exampleCAGEset,1 )
#' ))
#' 
#' cumulativeCTSSdistribution(exampleCAGEset, "consensusClusters") # Defaults in object do not fit
#' quantilePositions( object = exampleCAGEset, clusters = "consensusClusters"
#'                  , qLow = c(0.1, 0.2), qUp = c(0.8, 0.9))
#'                  
#' head(cbind(
#'   CAGEr:::consensusClustersQuantileLow(exampleCAGEset, 1),
#'   CAGEr:::consensusClustersQuantileUp (exampleCAGEset , 1)
#' ))
#' 
#' quantilePositions(exampleCAGEexp, "tagClusters",       qLow = c(0.1, 0.2), qUp = c(0.8, 0.9))
#' tagClustersGR(exampleCAGEexp)
#' quantilePositions(exampleCAGEexp, "consensusClusters", qLow = c(0.1, 0.2), qUp = c(0.8, 0.9))
#' 
#' @export

setGeneric( "quantilePositions"
          , function(object, clusters, qLow = 0.1, qUp = 0.9, useMulticore = FALSE, nrCores = NULL)
	standardGeneric("quantilePositions"))

#' @rdname quantilePositions

setMethod( "quantilePositions", "CAGEr"
         ,function (object, clusters, qLow, qUp, useMulticore, nrCores) {
	objName <- deparse(substitute(object))
	gr2tcq <- function(gr, q) {
	  tcq <- mcols(gr)[, paste("q", q, sep = "_"), drop = FALSE]
	  tcq <- data.frame(lapply(tcq, decode))
	  tcq <- tcq + start(gr)
	  cbind(cluster = names(gr), tcq)
	}
	message("\nGetting positions of quantiles within clusters...")
	if (clusters == "tagClusters") {
	  ctss.clusters <- lapply(sampleLabels(object), function(s) {
	    message("\t-> ", s)
	    .get.quant.pos( cumsums = CTSScumulativesTagClusters(object, s)
		                , clusters = tagClustersGR(object, s)
		                , q = c(qLow, qUp)
		                , useMulticore = useMulticore
		                , nrCores = nrCores)
	  })
	  names(ctss.clusters) <- sampleLabels(object)
	  ctss.clusters <- GRangesList(ctss.clusters)
		if(class(object) == "CAGEexp") {
		  tagClustersGR(object) <- ctss.clusters
		} else if (class(object) == "CAGEset") {
		  for (s in names(ctss.clusters)) {
	  	  tagClustersQuantileLow(object, s) <- gr2tcq(ctss.clusters[[s]], qLow)
        tagClustersQuantileUp (object, s) <- gr2tcq(ctss.clusters[[s]], qUp)
		  }
		} else stop("Unsupported CAGEr class.")
	} else if (clusters == "consensusClusters"){
	  cons.clusters.l <- lapply(sampleLabels(object), function(s) {
	    message("\t-> ", s)
	    .get.quant.pos( cumsums = CTSScumulativesCC(object, s)
		                , clusters = consensusClustersGR(object)
		                , q = c(qLow, qUp)
		                , useMulticore = useMulticore
		                , nrCores = nrCores)
	  })
	  names(cons.clusters.l) <- sampleLabels(object)
	  cons.clusters.l <- GRangesList(cons.clusters.l)
	  
		if(class(object) == "CAGEexp") {
      for (qName in paste("q", c(qLow, qUp), sep = "_")) {
        assays(consensusClustersSE(object))[[qName]] <-
          DataFrame(lapply(cons.clusters.l, function(gr) mcols(gr)[,qName]))}
		} else if (class(object) == "CAGEset") {
		  for (s in names(cons.clusters.l)) {
	  	  consensusClustersQuantileLow(object, s) <- gr2tcq(cons.clusters.l[[s]], qLow)
        consensusClustersQuantileUp (object, s) <- gr2tcq(cons.clusters.l[[s]], qUp)
		  }
		} else stop("Unsupported CAGEr class.")
	} else stop("'clusters' parameter must be one of the (\"tagClusters\", \"consensusClusters\")")
	assign(objName, object, envir = parent.frame())
	invisible(1)
})
