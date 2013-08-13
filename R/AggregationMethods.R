###################################################################
# Functions for aggregating tag clusters (TCs) across all samples
#

setGeneric(
name="aggregateTagClusters",
def=function(object, tpmThreshold = 5, excludeSignalBelowThreshold = TRUE, qLow = NULL, qUp = NULL, maxDist = 100){
	standardGeneric("aggregateTagClusters")
}
)

setMethod("aggregateTagClusters",
signature(object = "CAGEset"),
function (object, tpmThreshold = 5, excludeSignalBelowThreshold = TRUE, qLow = NULL, qUp = NULL, maxDist = 100){

	objName <- deparse(substitute(object))
	sample.labels = sampleLabels(object)	
	
	TC.list <- lapply(sample.labels, function(x) {tagClusters(object, sample = x)})
	names(TC.list) <- sample.labels
	if(length(qLow) > 0 & length(qUp) > 0){
		q.low.list <- object@tagClustersQuantileLow
		q.up.list <- object@tagClustersQuantileUp
		TC.list <- lapply(sample.labels, function(x) {
						  
							tc <- TC.list[[x]]
							q.low <- q.low.list[[x]]
							colnames(q.low) <- c("cluster", sub("q_", "q.low_", colnames(q.low)[-1], fixed = T))
							q.up <- q.up.list[[x]]
							colnames(q.up) <- c("cluster", sub("q_", "q.up_", colnames(q.up)[-1], fixed = T))
							tc <- merge(tc, q.low, by.x = "cluster", by.y = "cluster")
							tc <- merge(tc, q.up, by.x = "cluster", by.y = "cluster")
							return(tc)
						  
						  }
						)
		names(TC.list) <- sample.labels
		if(paste("q.low_", qLow, sep = "") %in% colnames(TC.list[[1]]) & paste("q.up_", qUp, sep = "") %in% colnames(TC.list[[1]])){
			consensus.clusters <- .make.consensus.clusters(TC.list = TC.list, start.coor = paste("q.low_", qLow, sep = ""), end.coor = paste("q.up_", qUp, sep = ""), plus.minus = round(maxDist/2), tpm.th = tpmThreshold)
		}else{
			stop("No data for given quantile positions! Run 'quantilePositions()' function for desired quantiles first, or omit 'qLow' and 'qUp' parameters to use start and end coordinates for aggregation instead!")
		}
	}else{
		consensus.clusters <- .make.consensus.clusters(TC.list = TC.list, start.coor = "start", end.coor = "end", plus.minus = round(maxDist/2), tpm.th = tpmThreshold)		
	}

	object@tagClustersInConsensusClusters <- consensus.clusters[,c("consensus.cluster", "cluster", "sample")]
	
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
		
		idx <- object@filteredCTSSidx
		ctss <- cbind(CTSScoordinates(object)[idx,], object@normalizedTpmMatrix[idx,,drop=F])

		tpm.list <- lapply(sample.labels, function(x) {ctss.s <- ctss[,c("chr", "pos", "strand", x)]; colnames(ctss.s)[4] <- "tagcount"; .getTotalTagCount(ctss.df = ctss.s, ctss.clusters = consensus.clusters, id.column = "consensus.cluster")})
		m <- matrix(unlist(tpm.list), nrow = nrow(consensus.clusters))
		dim.n <- list(consensus.cluster = as.character(consensus.clusters$consensus.cluster), sample = unname(sample.labels))
		dimnames(m) <- dim.n
		consensus.clusters$tpm <- rowSums(m) 
		
	}
	
	object@consensusClustersTpmMatrix <- m
	object@consensusClusters <- consensus.clusters
	
	assign(objName, object, envir = parent.frame())
	invisible(1)

}
)
