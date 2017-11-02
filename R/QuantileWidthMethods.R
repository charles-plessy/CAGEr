#' @include CAGEr.R

#' @name quantilePositions
#' 
#' @title Determine CTSS quantile positions within clusters
#' 
#' @description Calculates the positions of \dQuote{upper} and \dQuote{lower}
#' quantiles of CAGE signal along \emph{tag clusters} or
#' \emph{consensus clusters} in each sample of a  CAGEr object.
#' 
#' @param object A \code{\link{CAGEr}} object.
#' 
#' @param clusters Either \code{tagClusters} or \code{consensusClusters}.
#' 
#' @param qLow,qUp Which \dQuote{lower} or \dQuote{upper} quantiles should be
#'        calculated. Numeric vector of values in range \code{[0,1]}.
#' 
#' @param useMulticore Logical, should multicore be used.
#'        \code{useMulticore = TRUE} has only effect on Unix-like platforms.
#' 
#' @param nrCores Number of cores to use when \code{useMulticore = TRUE}.
#'        Default value \code{NULL} uses all detected cores.
#' 
#' @details The position of a \dQuote{lower} quantile \code{qLow} is defined as
#' a point that divides the cluster into two parts, so that the 5' part contains
#' \code{< qLow * 100\%} of the CAGE signal of that region (\code{>= qUp * 100\%}
#' for \dQuote{upper} quantiles).  Therefore, when \code{qLow <= qUp}, a lower
#' and and uppper quantile define a central region that contains
#' \code{>= (qUp - qLow) * 100\%} of the CAGE signal of that cluster.  The
#' width of that central part is refered to as "interquantile width", which is
#' a more robust measure of the promoter width than the total length of the
#' cluster.
#' 
#' @return When \code{clusters = "tagClusters"}, the slots \code{tagClustersQuantileLow}
#' and \code{tagClustersQuantileUp} of a provided \code{\link{CAGEset}} object will
#' be occupied with the positions of specified quantiles in all tag clusters for all
#' CAGE datasets. When \code{clusters = "consensusClusters"} the slots
#' \code{consensusClustersQuantileLow} and \code{consensusClustersQuantileUp} will be
#' occupied by the corresponding information for consensus clusters.
#' 
#' In \code{\link{CAGEexp}} objects, the positions of the quantiles are defined
#' reliatively to the start point of their cluster, for more efficient
#' \code{Rle} compression.  The quantile data for \emph{tag clusters} are stored
#' in the \code{TagClusters} objects directly.  The quantile data for
#' \emph{consensus clusters} are stored in \code{\link{integer}} matrices named
#' \dQuote{q_\emph{x}}, where \emph{x} represents the quantile (for instance,
#' \code{q_0.1}), and these matrices are \emph{assays} of the
#' \code{consensusClusters} \code{\link{RangedSummarizedExperiment}}.
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
          , function( object
                    , clusters = c("tagClusters", "consensusClusters")
                    , qLow = 0.1, qUp = 0.9
                    , useMulticore = FALSE, nrCores = NULL)
	standardGeneric("quantilePositions"))

#' @rdname quantilePositions

setMethod( "quantilePositions", "CAGEr"
         , function (object, clusters, qLow, qUp, useMulticore, nrCores) {
	objName <- deparse(substitute(object))
	gr2tcq <- function(gr, q) {
	  tcq <- mcols(gr)[, paste("q", q, sep = "_"), drop = FALSE]
	  tcq <- data.frame(lapply(tcq, decode))
	  tcq <- tcq + start(gr)
	  cbind(cluster = names(gr), tcq)
	}
	clusters <- match.arg(clusters)
	message("\nGetting positions of quantiles within clusters...")
	if (clusters == "tagClusters") {
	  ctss.clusters <- lapply(sampleList(object), function(s) {
	    message("\t-> ", s)
	    .get.quant.pos( cumsums = CTSScumulativesTagClusters(object, s)
		                , clusters = tagClustersGR(object, s)
		                , q = c(qLow, qUp)
		                , useMulticore = useMulticore
		                , nrCores = nrCores)
	  })
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
	  cons.clusters.l <- lapply(sampleList(object), function(s) {
	    message("\t-> ", s)
	    .get.quant.pos( cumsums = CTSScumulativesCC(object, s)
		                , clusters = consensusClustersGR(object)
		                , q = c(qLow, qUp)
		                , useMulticore = useMulticore
		                , nrCores = nrCores)
	  })
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
	}
	assign(objName, object, envir = parent.frame())
	invisible(1)
})
