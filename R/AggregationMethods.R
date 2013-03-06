###################################################################
# Functions for aggregating tag clusters (TCs) across all samples
#

setGeneric(
name="aggregateTagClusters",
def=function(object, tpmThreshold = 5, qLow = NULL, qUp = NULL, maxDist = 100){
	standardGeneric("aggregateTagClusters")
}
)

setMethod("aggregateTagClusters",
signature(object = "CAGEset"),
function (object, tpmThreshold, qLow, qUp, maxDist){

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
	
	m <- tapply(consensus.clusters$tpm, INDEX = list(consensus.cluster = consensus.clusters$consensus.cluster, sample = consensus.clusters$sample), FUN = sum)
	m[is.na(m)] <- 0 

	object@consensusClustersTpmMatrix <- m

# using data.table package
#	consensus.clusters <- data.table(consensus.clusters)
#	setkey(consensus.clusters, consensus.cluster)
#	consensus.clusters <- consensus.clusters[, list(chr[1], min(start), max(end), strand[1], sum(tpm)), by = consensus.cluster]
#	setnames(consensus.clusters, c("consensus.cluster", "chr", "start", "end", "strand", "tpm"))

#	consensus.clusters <- lapply(as.list(unique(consensus.clusters$consensus.cluster)), function(x) {
#									t <- subset(consensus.clusters, consensus.cluster == x, drop = F)
#									a <- data.frame(consensus.cluster = x, chr = t$chr[1], start = min(t$start), end = max(t$end), strand = t$strand[1], tpm = sum(t$tpm))
#									return(a)
#								 })
#	consensus.clusters <- do.call(rbind, consensus.clusters)
	
	
	
	clusters <- consensus.clusters[!duplicated(consensus.clusters$consensus.cluster), c("consensus.cluster", "chr", "strand")]
	start <- aggregate(consensus.clusters$start, by = list(consensus.clusters$consensus.cluster), FUN = min)
	end <- aggregate(consensus.clusters$end, by = list(consensus.clusters$consensus.cluster), FUN = max)[,2]
	tpm <- aggregate(consensus.clusters$tpm, by = list(consensus.clusters$consensus.cluster), FUN = sum)[,2]
	consensus.clusters.1 <- data.frame(consensus.cluster = as.integer(start[,1]), start = as.integer(start[,2]), end = as.integer(end), tpm = tpm)
	consensus.clusters <- merge(clusters, consensus.clusters.1)
	consensus.clusters <- consensus.clusters[,c("consensus.cluster", "chr", "start", "end", "strand", "tpm")]
	
	
	object@consensusClusters <- consensus.clusters
	
	assign(objName, object, envir = parent.frame())
	invisible(1)

}
)
