#' @include AllClasses.R Annotations.R CAGEexp.R CAGEr.R ClusteringMethods.R ClusteringFunctions.R CTSS.R Multicore.R

#' @name aggregateTagClusters
#' 
#' @title Aggregate TCs across all samples
#' 
#' @description Aggregates tag clusters (TCs) across all CAGE datasets within the CAGEr object
#' to create a referent set of consensus clusters.
#' 
#' @param object A \code{\link{CAGEr}} object
#' 
#' @param tpmThreshold Ignore tag clusters with normalized signal \code{< tpmThreshold} when
#'        constructing the consensus clusters.
#' 
#' @param excludeSignalBelowThreshold When \code{TRUE} the tag clusters with normalized signal
#'        \code{< tpmThreshold} will not contribute to the total CAGE signal of a consensus.
#'        cluster.  When set to \code{FALSE} all TCs that overlap consensus cluster will
#'        contribute to the total signal, regardless whether they pass the threshold for
#'        constructing the clusters or not.
#'        
#' @param qLow,qUp Set which "lower" (or "upper") quantile should be used as 5'
#'        (or 3') boundary of the tag cluster.  If \code{qLow = NULL} or \code{qUp = NULL},
#'        the start (or end) position of the TC is used.
#' 
#' @param maxDist Maximal length of the gap (in base-pairs) between two tag clusters for them to
#'        be part of the same consensus clusters.
#'        
#' @param useMulticore Logical, should multicore be used.  \code{useMulticore = TRUE} is
#'        supported only on Unix-like platforms.
#' 
#' @param nrCores Number of cores to use when \code{useMulticore = TRUE}.  Default value
#'        \code{NULL} uses all detected cores.
#' 
#' @details Since the tag clusters (TCs) returned by the \code{\link{clusterCTSS}} function
#' are constructed separately for every CAGE sample within the CAGEr object, they can differ
#' between samples in both their number, genomic coordinates, position  of dominant TSS and
#' overall signal.  To be able to compare all samples at the level of clusters of TSSs, TCs
#' from all CAGE datasets are aggregated into a single set of consensus clusters.
#' First, TCs with signal \code{>= tpmThreshold} from all CAGE datasets are selected, and their
#' 5' and 3' boundaries are determined based on provided \code{qLow} and \code{qUp} parameter
#' (or the start and end coordinates, if \code{qLow = NULL} and \code{qUp = NULL}.
#' Finally, the defined set of TCs from all CAGE datasets is reduced to a non-overlapping set
#' of consensus clusters by merging overlapping TCs and TCs \code{<= maxDist} base-pairs apart.
#' Consensus clusters represent a referent set of promoters that can be further used for
#' expression profiling or detecting "shifting" (differentially used) promoters between different
#' CAGE samples.
#' 
#' @return For \code{\link{CAGEset}} objects, the \code{consensusClusters} slot will be
#' populated with a data frame indicating the cluster name, chromosome, start and end
#' coordinates, the strand, and the normalised expression score of the cluster.  This
#' table is returned by the \code{\link{consensusClusters}} function.
#' 
#' For \code{\link{CAGEexp}} objects, the experiment \code{consensusClusters}
#' will be occupied by a \code{\link{RangedSummarizedExperiment}} containing the cluster
#' coodinates as row ranges, and their expression levels in the \code{counts} and \code{normalized}
#' assays.  These genomic ranges are returned by the \code{\link{consensusClustersGR}} function.
#' The CTSS ranges of the \code{tagCountMatrix} experiment will gain a
#' \code{cluster} column indicating which cluster they belong to.  Lastly, the number of
#' CTSS outside clusters will be documented in the \code{outOfClusters} column data.
#' This table is returned by the \code{\link{consensusClusters}} function.
#' 
#' @author Vanja Haberle
#' @author Charles Plessy
#' 
#' @family CAGEr object modifiers
#' @family CAGEr clusters functions
#' 
#' @importFrom data.table data.table setkeyv setnames
#' 
#' @examples
#' head(consensusClusters(exampleCAGEset))
#' aggregateTagClusters( exampleCAGEset, tpmThreshold = 50
#'                     , excludeSignalBelowThreshold = FALSE, maxDist = 100)
#' head(consensusClusters(exampleCAGEset))
#' 
#' aggregateTagClusters(object = exampleCAGEset, tpmThreshold = 50,
#'   excludeSignalBelowThreshold = FALSE, qLow = 0.1, qUp = 0.9, maxDist = 100)
#' head(consensusClusters(exampleCAGEset))
#' 
#' consensusClustersGR(exampleCAGEexp)
#' aggregateTagClusters( exampleCAGEexp, tpmThreshold = 50
#'                     , excludeSignalBelowThreshold = FALSE, maxDist = 100)
#' consensusClustersGR(exampleCAGEexp)
#' 
#' aggregateTagClusters( exampleCAGEexp, tpmThreshold = 50
#'                     , excludeSignalBelowThreshold = TRUE, maxDist = 100)
#' consensusClustersGR(exampleCAGEexp)
#' 
#' aggregateTagClusters( exampleCAGEexp, tpmThreshold = 50
#'                     , excludeSignalBelowThreshold = TRUE, maxDist = 100
#'                     , qLow = 0.1, qUp = 0.9)
#' consensusClustersGR(exampleCAGEexp)
#' 
#' @export

setGeneric( "aggregateTagClusters"
          , function( object
                    , tpmThreshold = 5, excludeSignalBelowThreshold = TRUE
                    , qLow = NULL, qUp = NULL
                    , maxDist = 100
                    , useMulticore = FALSE, nrCores = NULL)
          	  standardGeneric("aggregateTagClusters"))

#' @rdname aggregateTagClusters

setMethod( "aggregateTagClusters", "CAGEr"
         , function ( object, tpmThreshold, excludeSignalBelowThreshold
                    , qLow, qUp, maxDist, useMulticore, nrCores) {
  objname <- deparse(substitute(object))

  consensus.clusters <- .aggregateTagClustersGR( object, tpmThreshold = tpmThreshold
                                               , qLow = qLow, qUp = qUp, maxDist = maxDist)
  
  if (excludeSignalBelowThreshold) {
    filter <- .filterCtss( object
                         , threshold       = tpmThreshold
                         , nrPassThreshold = 1
                         , thresholdIsTpm  = TRUE)
  } else filter <- TRUE
	
	if (inherits(object, "CAGEset")) {
	  filteredSE <- CTSStagCountSE(object)[filter,]
	  .getTotalTagCountSample <- function(x) {
      ctss.s        <- .CTSS(rowRanges(filteredSE))
      score(ctss.s) <- assay(filteredSE, "normalizedTpmMatrix")[[x]]
      ctss.s        <- ctss.s[ctss.s$filteredCTSSidx]
      ctss.s        <- ctss.s[score(ctss.s) > 0]
  	  .getTotalTagCount(ctss = ctss.s, ctss.clusters = consensus.clusters)}
    tpm.list <- bplapply( sampleList(object), .getTotalTagCountSample
                        , BPPARAM = CAGEr_Multicore(useMulticore, nrCores))
		m <- as.matrix(data.frame(tpm.list))
		rownames(m) <- 1:nrow(m)
		consensus.clusters$tpm <- rowSums(m)
	  object@consensusClustersTpmMatrix <- m
	  consensusClusters(object) <- CCgranges2dataframe(consensus.clusters)
	}
	
	if (inherits(object, "CAGEexp")) {
    CTSScoordinatesGR(object)$cluster <-
      ranges2names(CTSScoordinatesGR(object), consensus.clusters)
    se <- CTSStagCountSE(object)[filter & decode(filteredCTSSidx(object)), ]
    consensusClustersSE(object) <- .CCtoSE(se, consensus.clusters)
    score(consensusClustersGR(object)) <- rowSums(assays(consensusClustersSE(object))[["normalized"]])
    object$outOfClusters <- librarySizes(object) - colSums(assay(consensusClustersSE(object)))
	}
	
	assign(objname, object, envir = parent.frame())
	invisible(1)
})


setGeneric( ".aggregateTagClustersGR"
          , function( object, tpmThreshold = 5
                    , qLow = NULL, qUp = NULL, maxDist = 100)
          	  standardGeneric(".aggregateTagClustersGR"))

setMethod( ".aggregateTagClustersGR", "CAGEr"
         , function ( object, tpmThreshold
                    , qLow, qUp, maxDist) {
  if (all( !is.null(qLow), !is.null(qUp))) {
    TC.list <- tagClustersGR(object, returnInterquantileWidth = TRUE,  qLow = qLow, qUp = qUp)
    TC.list <- endoapply(TC.list, function(x) {
      start(x) <- mcols(x)[[paste0("q_", qLow)]] + start(x)
      end(x)   <- mcols(x)[[paste0("q_", qUp)]]  + start(x)
      x})
  } else {
    TC.list <- tagClustersGR(object)
  }
  consensus.clusters <- .make.consensus.clusters( TC.list = TC.list
                                                , plus.minus = round(maxDist/2)
                                                , tpm.th = tpmThreshold)
  consensus.clusters <- .clusterAggregateAndSum(consensus.clusters, "consensus.cluster")
  consensus.clusters <- GRanges(consensus.clusters)
  names(consensus.clusters) <- as.character(consensus.clusters)
  .ConsensusClusters(consensus.clusters)
})


setGeneric( ".CCtoSE" , function(se, consensus.clusters, tpmThreshold = 1)
          	  standardGeneric(".CCtoSE"))

setMethod( ".CCtoSE"
         , c(se = "RangedSummarizedExperiment")
         , function(se, consensus.clusters, tpmThreshold = 1) {
    if (is.null(assays(se)[["normalizedTpmMatrix"]]))
      stop("Needs normalised data; run ", sQuote("normalizeTagCount()"), " first.")
    if (is.null(rowRanges(se)$cluster))
      rowRanges(se)$cluster <- ranges2names(rowRanges(se), consensus.clusters)
    
    if (tpmThreshold > 0)
      se <- se[rowSums(DelayedArray(assays(se)[["normalizedTpmMatrix"]])) > tpmThreshold,]
    
    .rowsumAsMatrix <- function(DF, names) {
      rs <- rowsum(as.matrix(DelayedArray(DF)), as.factor(names))
      if (rownames(rs)[1] == "") # If some CTSS were not in clusters
        rs <- rs[-1,]
      rs
      }
    
    counts <- .rowsumAsMatrix(assays(se)[["counts"]], rowRanges(se)$cluster)
    norm   <- .rowsumAsMatrix(assays(se)[["normalizedTpmMatrix"]], rowRanges(se)$cluster)

	  SummarizedExperiment( rowRanges = consensus.clusters[rownames(counts)]
	                      , assays    = SimpleList( counts     = counts
	                                              , normalized = norm))
})

#' @name CustomConsensusClusters
#' 
#' @title Expression levels of consensus cluster
#' 
#' @description Intersects custom consensus clusters with the CTSS data in a 
#' [`CAGEexp`] object, and stores the result as a expression matrices
#' (raw and normalised tag counts).
#' 
#' @param object A `CAGEexp` object
#' 
#' @param clusters Consensus clusters in [`GRanges`] format.
#' 
#' @param threshold,nrPassThreshold Only CTSSs with signal `>= threshold` in
#'        `>= nrPassThreshold` experiments will be used for clustering and will
#'        contribute towards total signal of the cluster.
#' 
#' @param thresholdIsTpm Logical, is threshold raw tag count value (FALSE) or
#'        normalized signal (TRUE).
#'        
#' @details Consensus clusters must not overlap, so that a single base of the
#' genome can only be attributed to a single cluster.  This is enforced by the
#' [.ConsensusCluster()] constructor.
#' 
#' @return stores the result as a new [`RangedSummarizedExperiment`] in the
#' `experiment` slot of the object.  The assays of the new experiment are called
#' `counts` and `normalized`.  An `outOfClusters` column is added
#' to the sample metadata to reflect the number of molecules that do not have
#' their TSS in a consensus cluster.
#' 
#' @author Charles Plessy
#' 
#' @family CAGEr object modifiers
#' @family CAGEr clusters functions
#'  
#' @examples 
#' 
#' cc <- consensusClustersGR(exampleCAGEexp)
#' CustomConsensusClusters(exampleCAGEexp, cc)
#' 
#' @export

setGeneric( "CustomConsensusClusters"
          , function( object
                    , clusters
                    , threshold       = 0
                    , nrPassThreshold = 1
                    , thresholdIsTpm  = TRUE)
          	  standardGeneric("CustomConsensusClusters"))

#' @rdname CustomConsensusClusters

setMethod( "CustomConsensusClusters", c("CAGEexp", "GRanges")
         , function (object, clusters
                    , threshold, nrPassThreshold, thresholdIsTpm  = TRUE) {
  objname <- deparse(substitute(object))
  
  clusters <- .ConsensusClusters(clusters)
  
  filter <- .filterCtss( object
                       , threshold       = threshold
                       , nrPassThreshold = nrPassThreshold
                       , thresholdIsTpm  = thresholdIsTpm)
  
  CTSScoordinatesGR(object)$cluster <- ranges2names(CTSScoordinatesGR(object), clusters)
  
  consensusClustersSE(object) <- .CCtoSE( CTSStagCountSE(object)[filter, ]
                                        , clusters)
  score(consensusClustersGR(object)) <- rowSums(assays(consensusClustersSE(object))[["normalized"]])
  object$outOfClusters <- librarySizes(object) - colSums(assay(consensusClustersSE(object)))
  
  assign(objname, object, envir = parent.frame())
	invisible(object)
})