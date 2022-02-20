#' @include CAGEr.R

#' @name quantilePositions
#' 
#' @title Determine CTSS quantile positions within clusters
#' 
#' @description Calculates the positions of \dQuote{upper} and \dQuote{lower}
#' quantiles of CAGE signal along _tag clusters_ or _consensus clusters_ in each
#' sample of a [`CAGEexp`] object.
#' 
#' @param object A CAGEexp object.
#' 
#' @param clusters Either `tagClusters` or `consensusClusters`.
#' 
#' @param qLow,qUp Which \dQuote{lower} or \dQuote{upper} quantiles should be
#'        calculated. Numeric vector of values in range `[0,1]`.
#' 
#' @param useMulticore Logical, should multicore be used.  `useMulticore = TRUE`
#'        has only effect on Unix-like platforms.
#' 
#' @param nrCores Number of cores to use when `useMulticore = TRUE`.  Default
#'        value `NULL` uses all detected cores.
#' 
#' @details From the 5' end the position, the position of a quantile _q_ is
#' determined as the first base in which of the cumulative expression is higher
#' or equal to _q%_ of the total CAGE signal of that cluster.  Promoter
#' _interquantile width_ is defined as the distance (in base pairs) between a
#' \dQuote{lower} and an \dQuote{upper} quantile position.
#' 
#' @return Returns the objects, in which the positions of the quantiles are
#' defined relatively to the start point of their cluster, for more efficient
#' `Rle` compression.  The quantile data for _tag clusters_ are stored in the
#' `TagClusters` objects directly.  The quantile data for `consensus clusters`
#' are stored in [`integer`] matrices named \dQuote{q_\emph{x}}, where \emph{x}
#' represents the quantile (for instance, `q_0.1`), and these matrices are
#' _assays_ of the `consensusClusters` [`RangedSummarizedExperiment`].
#' 
#' @author Vanja Haberle
#' @author Charles Plessy
#' 
#' @family CAGEr object modifiers
#' @family CAGEr clusters functions
#' 
#' @examples 
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

setMethod( "quantilePositions", "CAGEexp"
         , function (object, clusters, qLow, qUp, useMulticore, nrCores) {
  clusters <- match.arg(clusters)
  message("\nGetting positions of quantiles within clusters...")
  if (clusters == "tagClusters") {
    ctss.clusters <- bplapply(sampleList(object), function(s) {
        message("\t-> ", s)
        .get.quant.pos( cum.sums  = CTSScumulativesTagClusters(object, s)
                      , clusters = tagClustersGR(object, s)
                      , q        = c(qLow, qUp))
      },
      BPPARAM = CAGEr_Multicore(useMulticore, nrCores))
    ctss.clusters <- GRangesList(ctss.clusters)
    tagClustersGR(object) <- ctss.clusters
  } else if (clusters == "consensusClusters"){
    cons.clusters.l <- bplapply(sampleList(object), function(s) {
        message("\t-> ", s)
        .get.quant.pos( cum.sums  = CTSScumulativesCC(object, s)
                      , clusters = consensusClustersGR(object)
                      , q        = c(qLow, qUp))
      },
      BPPARAM = CAGEr_Multicore(useMulticore, nrCores))
    for (qName in paste("q", c(qLow, qUp), sep = "_")) {
      assays(consensusClustersSE(object), withDimnames=FALSE)[[qName]] <-
        DataFrame(lapply(cons.clusters.l, function(gr) mcols(gr)[,qName]))
    }
  }
  object
})

#' Get quantile positions
#' 
#' Private function that calculates position of quantiles for CTSS clusters
#' based on distribution of tags within the clusters.
#' 
#' @param cum.sums Named list of vectors containing cumulative sum for each
#'        cluster (returned by the `CTSScumulativesTagClusters` or
#'        `CTSScumulativesCC` function).
#' @param clusters [`TagClusters`] or [`ConsensusClusters`] object representing
#'        tag clusters or consensus clusters.
#' @param q desired quantiles - single value or a vector of values.
#' 
#' @return Returns the `clusters` object with one more metadata column per value
#' in `q`, containing `Rle` integers giving the relative distance of the
#' quantile boundaries to the start position.
#' 
#' @examples 
#' cum.sums  <- RleList(`1` = Rle(1), `2` = cumsum(Rle(c(1, 1, 1, 2, 4, 0, 1, 1))))
#' clusters <- GRanges(c("chr1:100-101", "chr1:120-127"))
#' CAGEr:::.get.quant.pos(cum.sums, clusters, c(.2, .8))

#' @rdname QuantileWidthFunctions

.get.quant.pos <- function(cum.sums, clusters, q) {
  # Vectorized function calculating one quantile position
  # for each element of a list of cumulative sums.
  getQuantilepos <- Vectorize(vectorize.args = "cum.sum", function(q, cum.sum) {
    cum.sum <- decode(cum.sum) # A microbenchmark showed it it 3 times faster when applying decode() now
    c.max <- tail(cum.sum,1) # Max is last element since x is a cumulative sums.
    treshold <- c.max * q
    which.max(cum.sum >= treshold)
  })
  # Calculate quantile positions for each quantile.
  cluster.q <- lapply(q, getQuantilepos, cum.sums)
  names(cluster.q) = paste('q_', q, sep = '')
  # Add one metadata column per quantile to the cluster object and return it.
  mcols(clusters)[, names(cluster.q)] <- DataFrame(lapply(cluster.q, Rle))
  clusters
}