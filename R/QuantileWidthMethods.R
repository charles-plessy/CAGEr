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
#' @details From the 5' end the position, the position of a quantile \emph{q} is
#' detemined as the first base in which of the cumulative expression is higher
#' or equal to \emph{q%} of the total CAGE signal of that cluster.  Promoter
#' \emph{interquantile width} is defined as the distance (in base pairs)
#' between a \dQuote{lower} and an \dQuote{upper} quantile position.
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
	clusters <- match.arg(clusters)
	message("\nGetting positions of quantiles within clusters...")
	if (clusters == "tagClusters") {
    ctss.clusters <- bplapply(sampleList(object),
      function(s) {
        message("\t-> ", s)
        .get.quant.pos( cumsums  = CTSScumulativesTagClusters(object, s)
                      , clusters = tagClustersGR(object, s)
                      , q        = c(qLow, qUp))
      },
      BPPARAM = CAGEr_Multicore(useMulticore, nrCores))
	  ctss.clusters <- GRangesList(ctss.clusters)
		if(inherits(object, "CAGEexp")) {
		  tagClustersGR(object) <- ctss.clusters
		} else if (inherits(object,  "CAGEset")) {
		  for (s in names(ctss.clusters)) {
	  	  tagClustersQuantileLow(object, s) <- tagClustersQuantile(ctss.clusters[[s]], q = qLow)
        tagClustersQuantileUp (object, s) <- tagClustersQuantile(ctss.clusters[[s]], q = qUp)
		  }
		} else stop("Unsupported CAGEr class.")
	} else if (clusters == "consensusClusters"){
    cons.clusters.l <- bplapply(sampleList(object),
      function(s) {
        message("\t-> ", s)
        .get.quant.pos( cumsums  = CTSScumulativesCC(object, s)
                      , clusters = consensusClustersGR(object)
                      , q        = c(qLow, qUp))
      },
      BPPARAM = CAGEr_Multicore(useMulticore, nrCores))
    
		if(inherits(object, "CAGEexp")) {
      for (qName in paste("q", c(qLow, qUp), sep = "_")) {
        assays(consensusClustersSE(object), withDimnames=FALSE)[[qName]] <-
          DataFrame(lapply(cons.clusters.l, function(gr) mcols(gr)[,qName]))}
		} else if (inherits(object, "CAGEset")) {
		  for (s in names(cons.clusters.l)) {
	  	  consensusClustersQuantileLow(object, s) <- .getClustersQuantile(cons.clusters.l[[s]], q = qLow)
        consensusClustersQuantileUp (object, s) <- .getClustersQuantile(cons.clusters.l[[s]], q = qUp)
		  }
		} else stop("Unsupported CAGEr class.")
	}
	assign(objName, object, envir = parent.frame())
	invisible(1)
})
