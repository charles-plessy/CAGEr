#' @include AllClasses.R Annotations.R CAGEexp.R CAGEr.R ClusteringMethods.R ClusteringFunctions.R CTSS.R Multicore.R

#' @name aggregateTagClusters
#' 
#' @title Aggregate TCs across all samples
#' 
#' @description Aggregates tag clusters (TCs) across all CAGE datasets within
#' the CAGEr object to create a referent set of consensus clusters.
#' 
#' @param object A [`CAGEr`] object
#' 
#' @param tpmThreshold Ignore tag clusters with normalized signal `< tpmThreshold`
#'        when constructing the consensus clusters.
#' 
#' @param excludeSignalBelowThreshold When `TRUE` the tag clusters with
#'        normalized signal `< tpmThreshold` will not contribute to the total
#'        CAGE signal of a consensus cluster.  When set to `FALSE` all TCs that
#'        overlap consensus clusters will contribute to the total signal,
#'        regardless whether they pass the threshold for constructing the
#'        clusters or not.
#'        
#' @param qLow,qUp Set which "lower" (or "upper") quantile should be used as 5'
#'        (or 3') boundary of the tag cluster.  If `NULL` the start (for `qLow`)
#'        or end (for `qUp`) position of the TC is used.
#' 
#' @param maxDist Maximal length of the gap (in base-pairs) between two tag
#'        clusters for them to be part of the same consensus clusters.
#'        
#' @param useMulticore Logical, should multicore be used (supported only on
#'        Unix-like platforms).
#' 
#' @param nrCores Number of cores to use when `useMulticore = TRUE`.  Default
#'        (`NULL`) uses all detected cores.
#' 
#' @details Since the tag clusters (TCs) returned by the [`clusterCTSS`]
#' function are constructed separately for every CAGE sample within the CAGEr
#' object, they can differ between samples in both their number, genomic
#' coordinates, position of dominant TSS and overall signal.  To be able to
#' compare all samples at the level of clusters of TSSs, TCs from all CAGE
#' datasets are aggregated into a single set of consensus clusters.  First, TCs
#' with signal `>= tpmThreshold` from all CAGE datasets are selected, and their
#' 5' and 3' boundaries are determined based on provided `qLow` and `qUp`
#' parameter (or the start and end coordinates, if they are set to `NULL`).
#' Finally, the defined set of TCs from all CAGE datasets is reduced to a
#' non-overlapping set of consensus clusters by merging overlapping TCs and TCs
#' `<= maxDist` base-pairs apart.  Consensus clusters represent a referent set
#' of promoters that can be further used for expression profiling or detecting
#'  "shifting" (differentially used) promoters between different CAGE samples.
#' 
#' @return Returns the object in which the _experiment_ `consensusClusters` will
#' be occupied by a [`RangedSummarizedExperiment`] containing the cluster
#' coordinates as row ranges, and their expression levels in the `counts` and
#' `normalized` assays.  These genomic ranges are returned by the
#' [`consensusClustersGR`] function and the whole object can be accessed with
#' the [`consensusClustersSE`] function.  The CTSS ranges of the
#' `tagCountMatrix` _experiment_ will gain a `cluster` column indicating which
#' cluster they belong to.  Lastly, the number of CTSS outside clusters will be
#' documented in the `outOfClusters` column data.
#'  
#' @author Vanja Haberle
#' @author Charles Plessy
#' 
#' @family CAGEr object modifiers
#' @family CAGEr clusters functions
#' 
#' @importFrom IRanges reduce
#' @importFrom GenomicRanges granges
#' @importFrom S4Vectors endoapply mcols
#' 
#' @examples
#' 
#' consensusClustersGR(exampleCAGEexp)
#' ce <- aggregateTagClusters( exampleCAGEexp, tpmThreshold = 50
#'                           , excludeSignalBelowThreshold = FALSE, maxDist = 100)
#' consensusClustersGR(ce)
#' 
#' ce <- aggregateTagClusters( exampleCAGEexp, tpmThreshold = 50
#'                           , excludeSignalBelowThreshold = TRUE, maxDist = 100)
#' consensusClustersGR(ce)
#' 
#' ce <- aggregateTagClusters( exampleCAGEexp, tpmThreshold = 50
#'                           , excludeSignalBelowThreshold = TRUE, maxDist = 100
#'                           , qLow = 0.1, qUp = 0.9)
#' consensusClustersGR(ce)
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
	
    CTSScoordinatesGR(object)$cluster <-
      ranges2names(CTSScoordinatesGR(object), consensus.clusters)
    se <- CTSStagCountSE(object)[filter & decode(filteredCTSSidx(object)), ]
    consensusClustersSE(object) <- .CCtoSE(se, consensus.clusters, tpmThreshold = tpmThreshold)
    score(consensusClustersGR(object)) <- rowSums(assays(consensusClustersSE(object))[["normalized"]])
    object$outOfClusters <- librarySizes(object) - colSums(assay(consensusClustersSE(object)))
    object
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
      end(x)   <- mcols(x)[[paste0("q_", qUp) ]] + start(x) - 1L
      start(x) <- mcols(x)[[paste0("q_", qLow)]] + start(x) - 1L
      x})
  } else {
    TC.list <- tagClustersGR(object)
  }

  # Filter out TCs with too low score.
  gr.list <- endoapply(TC.list, function (gr) gr <- gr[score(gr) >= tpmThreshold])
  
  # Aggregate clusters by expanding and merging TCs from all samples.
  clusters.gr <- unlist(gr.list)
  clusters.gr <- reduce(clusters.gr, min.gapwidth = maxDist)
  
  # CTSS with score that is sum od all samples
  ctss <- CTSScoordinatesGR(object)
  score(ctss) <- rowSums(CTSSnormalizedTpmDF(object) |> DelayedArray::DelayedArray() )
  
  # See `benchmarks/dominant_ctss.md`.
  o <- findOverlaps(clusters.gr, ctss)
  
  rl <- rle(queryHits(o))$length
  cluster_start_idx <- cumsum(c(1, head(rl, -1))) # Where each run starts
  grouped_scores <- extractList(score(ctss), o)
  # grouped_pos    <- extractList(pos(ctss), o)
  
  find.dominant.idx <- function (x) {
    # which.max is breaking ties by taking the last, but this will give slightly
    # different biases on plus an minus strands.
    w <- which(x == max(x))
    w[ceiling(length(w)/2)]
  }
  local_max_idx <- sapply(grouped_scores, find.dominant.idx) -1  # Start at zero
  global_max_ids <- cluster_start_idx + local_max_idx
  # start(clusters.gr) <- min(grouped_pos)
  # end  (clusters.gr) <- max(grouped_pos)
  score(clusters.gr) <- sum(grouped_scores)
  clusters.gr$dominant_ctss <- granges(ctss)[subjectHits(o)][global_max_ids]
  clusters.gr$tpm.dominant_ctss <- score(ctss)[subjectHits(o)][global_max_ids]
  clusters.gr

  names(clusters.gr) <- as.character(clusters.gr)
  .ConsensusClusters(clusters.gr)
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
      # First, remove CTSS that do not match clusters
      DF <- DF[names != "",]
      names <- names[names != ""]
      rowsum(as.matrix(DelayedArray(DF)), as.factor(names), reorder = FALSE)
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
#' [`.ConsensusClusters`] constructor.
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

  object
})
