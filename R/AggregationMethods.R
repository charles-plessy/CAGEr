#' @include AllClasses.R CAGEexp.R CAGEr.R ClusteringFunctions.R

#' @name aggregateTagClusters
#' 
#' @title Aggregate TCs across all samples
#' 
#' @description Aggregates tag clusters (TCs) across all CAGE datasets within the CAGEr object
#' to create a referent set of consensus clusters.
#' 
#' @param object A \code{\link{CAGEset}} object
#' 
#' @param tpmThreshold Only tag clusteres with normalized signal \code{>= tpmThreshold} will be
#'        used to construct consensus clusters.
#' 
#' @param excludeSignalBelowThreshold 	When \code{TRUE} only tag clusters with normalized signal
#'        \code{>= tpmThreshold} will contribute to the total CAGE signal of a consensus cluster,
#'        \emph{i.e.} only the TCs that are used to construct consensus cluster.  When set to
#'        \code{FALSE} all TCs that overlap consensus cluster will contribute to the total signal
#'        (regardless whether they pass the threshold or not), however only the TCs above the
#'        threshold will be used to define consensus cluster boundaries. Thus, it that case the TCs
#'        above the threshold are first used to construct consensus clusters and define their
#'        boundaries, but then CAGE signal from all TCs that fall within those boundaries is used
#'        to calculate total signal of a particular consensus cluster.
#'        
#' @param qLow,qUp Position of which "lower" (resp. "upper") quantile should be used as 5'
#'        (resp. 3') boundary of the tag cluster.  If \code{qLow = NULL} start position of the TC
#'        is used.  If \code{qUp = NULL} end position of the TC is used. \code{qUp} has to be
#'        \code{>= qLow}. See Details.
#' 
#' @param maxDist Maximal length of the gap (in base-pairs) between two tag clusters for them to
#'        be part of the same consensus clusters.  See Details.
#' 
#' @details Tag clusters (TCs) returned by \code{\link{clusterCTSS}} function are constructed for
#' every CAGE dataset within CAGEr object separatelly, based on the CAGE signal in that sample.
#' Thus, TCs from two CAGE datasets can differ both in their number, genomic coordinates, position
#' of dominant TSS and overall signal. To be able to compare all samples at the level of clusters
#' of TSSs, TCs from all CAGE datasets are aggregated into a single set of consensus clusters.
#' First, TCs with signal \code{>= tpmThreshold} from all CAGE datasets are selected, and their 5'
#' and 3' boundaries are determined based on provided \code{qLow} and \code{qUp} parameters. If
#' \code{qLow = NULL} and \code{qUp = NULL} the start and end coordinates, \emph{i.e.} the full
#' span of the TC is used, otherwise the positions of \code{qLow} and \code{qUp} quantiles are
#' used as 5' and 3' boundary, respectively.  Finally, the defined set of TCs from all CAGE
#' datasets is reduced to a non-overlapping set of consensus clusters by merging overlapping TCs
#' and TCs \code{<= maxDist} base-pairs apart.  Consensus clusters represent a referent set of
#' promoters that can be further used for expression profiling or detecting "shifting"
#' (differentially used) promoters between different CAGE samples.
#' 
#' @return For \code{\link{CAGEexp}} objects, the slots \code{consensusClusters},
#' \code{tagClustersInConsensusClusters} and \code{consensusClustersTpmMatrix} will be occupied
#' by the genomic coordinates of consensus clusters, information on containing TCs and the total
#' CAGE signal across all CAGE datasets, respectively.
#' 
#' @author Vanja Haberle
#' 
#' @family CAGEr object modifiers
#' @family CAGEr clusters functions
#' 
#' @importFrom data.table data.table setkey setnames
#' 
#' @examples
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' head(consensusClusters(exampleCAGEset))
#' aggregateTagClusters(object = exampleCAGEset, tpmThreshold = 50,
#'   excludeSignalBelowThreshold = FALSE, qLow = 0.1, qUp = 0.9, maxDist = 100)
#' head(consensusClusters(exampleCAGEset))
#' 
#' @export

setGeneric(
name="aggregateTagClusters",
def=function(object, tpmThreshold = 5, excludeSignalBelowThreshold = TRUE, qLow = NULL, qUp = NULL, maxDist = 100){
	standardGeneric("aggregateTagClusters")
}
)

setMethod("aggregateTagClusters",
signature(object = "CAGEr"),
function (object, tpmThreshold, excludeSignalBelowThreshold, qLow, qUp, maxDist){
	objName <- deparse(substitute(object))

	TC.list <- getTagCluster(object)
	if(length(qLow) > 0 & length(qUp) > 0){
		TC.list <- lapply(sampleLabels(object), function(x) {
							tc <- TC.list[[x]]
							q.low <- tagClustersQuantileLow(object, x)
							colnames(q.low) <- c("cluster", sub("q_", "q.low_", colnames(q.low)[-1], fixed = T))
							q.up <- tagClustersQuantileUp(object, x)
							colnames(q.up) <- c("cluster", sub("q_", "q.up_", colnames(q.up)[-1], fixed = T))
							tc <- merge(tc, q.low, by.x = "cluster", by.y = "cluster")
							tc <- merge(tc, q.up, by.x = "cluster", by.y = "cluster")
							return(tc)
						  })
		names(TC.list) <- sampleLabels(object)
		if(paste("q.low_", qLow, sep = "") %in% colnames(TC.list[[1]]) & paste("q.up_", qUp, sep = "") %in% colnames(TC.list[[1]])){
			consensus.clusters <- .make.consensus.clusters(TC.list = TC.list, start.coor = paste("q.low_", qLow, sep = ""), end.coor = paste("q.up_", qUp, sep = ""), plus.minus = round(maxDist/2), tpm.th = tpmThreshold)
		}else{
			stop("No data for given quantile positions! Run 'quantilePositions()' function for desired quantiles first, or omit 'qLow' and 'qUp' parameters to use start and end coordinates for aggregation instead!")
		}
	}else{
		consensus.clusters <- .make.consensus.clusters(TC.list = TC.list, start.coor = "start", end.coor = "end", plus.minus = round(maxDist/2), tpm.th = tpmThreshold)		
	}

	TCsInConsensusClusters(object) <- consensus.clusters[,c("consensus.cluster", "cluster", "sample")]
	
	if(excludeSignalBelowThreshold){
		m <- tapply(consensus.clusters$tpm, INDEX = list(consensus.cluster = consensus.clusters$consensus.cluster, sample = consensus.clusters$sample), FUN = sum)
		m[is.na(m)] <- 0
	}
	
	consensus.clusters <- data.table(consensus.clusters)
	setkey(consensus.clusters, consensus.cluster)
	consensus.clusters <- consensus.clusters[, list(chr[1], min(start), max(end), strand[1], sum(tpm)), by = consensus.cluster]
	setnames(consensus.clusters, c("consensus.cluster", "chr", "start", "end", "strand", "tpm"))
	setkey(consensus.clusters, consensus.cluster)
	consensus.clusters <- as.data.frame(consensus.clusters)
		
	if(!excludeSignalBelowThreshold){
		idx <- filteredCTSSidx(object)
		ctss <- cbind(CTSScoordinates(object)[idx,], CTSSnormalizedTpm(object)[idx,,drop=F])

		tpm.list <- lapply(sampleLabels(object), function(x) {ctss.s <- ctss[,c("chr", "pos", "strand", x)]; colnames(ctss.s)[4] <- "tagcount"; .getTotalTagCount(ctss.df = ctss.s, ctss.clusters = consensus.clusters, id.column = "consensus.cluster")})
		m <- matrix(unlist(tpm.list), nrow = nrow(consensus.clusters))
		dim.n <- list(consensus.cluster = as.character(consensus.clusters$consensus.cluster), sample = unname(sampleLabels(object)))
		dimnames(m) <- dim.n
		consensus.clusters$tpm <- rowSums(m) 
	}
	
	if (class(object) == "CAGEset") {
	  object@consensusClustersTpmMatrix <- m
	  consensusClusters(object) <- consensus.clusters
	}
	
	if (class(object) == "CAGEexp") {
	  consensusClustersSE(object) <-
	    SummarizedExperiment( rowRanges = CCdataframe2granges(consensus.clusters)
	                        , assay     = SimpleList(counts = m))
	}
	
	assign(objName, object, envir = parent.frame())
	invisible(1)
}
)
