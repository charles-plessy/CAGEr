#' @include AllClasses.R CAGEexp.R CAGEr.R ClusteringMethods.R ClusteringFunctions.R CTSS.R Multicore.R

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
#' @importFrom data.table data.table setkey setnames
#' 
#' @examples
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' head(consensusClusters(exampleCAGEset))
#' aggregateTagClusters( exampleCAGEset, tpmThreshold = 50
#'                     , excludeSignalBelowThreshold = FALSE, maxDist = 100)
#' head(consensusClusters(exampleCAGEset))
#' aggregateTagClusters(object = exampleCAGEset, tpmThreshold = 50,
#'   excludeSignalBelowThreshold = FALSE, qLow = 0.1, qUp = 0.9, maxDist = 100)
#' head(consensusClusters(exampleCAGEset))
#' 
#' ce <- readRDS(system.file(package = "CAGEr", "extdata/CAGEexp.rds"))
#' aggregateTagClusters(ce, tpmThreshold = 50, excludeSignalBelowThreshold = FALSE, maxDist = 100)
#' consensusClustersGR(ce)
#' aggregateTagClusters(ce, tpmThreshold = 50, excludeSignalBelowThreshold = TRUE, maxDist = 100)
#' consensusClustersGR(ce)
#' aggregateTagClusters( ce, tpmThreshold = 50, excludeSignalBelowThreshold = TRUE, maxDist = 100
#'                     , qLow = 0.1, qUp = 0.9)
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
	objName <- deparse(substitute(object))

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
  consensus.clusters <- .ConsensusClusters(consensus.clusters)
	if (excludeSignalBelowThreshold) {
		m <- tapply(score(consensus.clusters), INDEX = list(consensus.cluster = consensus.clusters$consensus.cluster, sample = consensus.clusters$sample), FUN = sum)
		m[is.na(m)] <- 0
		m <- m[, sampleLabels(object)]
	}

	consensus.clusters <- .clusterAggregateAndSum(consensus.clusters, "consensus.cluster")
	
  if (!excludeSignalBelowThreshold) {
    useMulticore <- .checkMulticore(useMulticore)
    .getTotalTagCountSample <- function(x) {
		  ctss.s <- CTSSnormalizedTpmGR(object, x)
		  ctss.s <- ctss.s[ctss.s$filteredCTSSidx]
		  .getTotalTagCount(ctss = ctss.s, ctss.clusters = consensus.clusters)}
    if(useMulticore == TRUE){
      if(is.null(nrCores))  nrCores <- detectCores()
      tpm.list <- mclapply(sampleLabels(object), .getTotalTagCountSample, mc.cores = nrCores)  
    }else{
      tpm.list <- lapply(sampleLabels(object), .getTotalTagCountSample)
    }
		names(tpm.list) <- sampleLabels(object)
		m <- as.matrix(data.frame(tpm.list))
		consensus.clusters$tpm <- rowSums(m) 
	}
	
	if (class(object) == "CAGEset") {
	  object@consensusClustersTpmMatrix <- m
	  consensusClusters(object) <- CCgranges2dataframe(consensus.clusters)
	}
	
	if (class(object) == "CAGEexp") {
    # See ranges2genes for similar code.
    names(consensus.clusters) <- as.character(consensus.clusters)
    cnames <- findOverlaps(CTSScoordinatesGR(object), consensus.clusters)
    cnames <- as(cnames, "List")
    cnames <- extractList(names(consensus.clusters), cnames)
    cnames <- unique(cnames)
    cnames <- unstrsplit(cnames, ";")
    CTSScoordinatesGR(object)$cluster <- Rle(cnames)
    counts <- rowsum(CTSStagCountDf(object), cnames)
    if (rownames(counts[1,]) == "") {  # If some CTSS were not in clusters
      object$outOfClusters <- unlist(counts[1,])
      counts <- counts[-1,]
    } else {
      object$outOfClusters <- 0
    }
    counts <- counts[names(consensus.clusters),]
    rownames(m)[consensus.clusters$consensus.cluster] <- names(consensus.clusters)
    consensusClustersSE(object) <-
	    SummarizedExperiment( rowRanges = consensus.clusters
	                        , assays    = SimpleList( counts = as.matrix(counts)
	                                                , normalized = m))
	}
	
	assign(objName, object, envir = parent.frame())
	invisible(1)
})