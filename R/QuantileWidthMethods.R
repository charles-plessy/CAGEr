setGeneric(
name="quantilePositions",
def=function(object, clusters, qLow = 0.1, qUp = 0.9, useMulticore = FALSE, nrCores = NULL){
	standardGeneric("quantilePositions")
}
)

setMethod("quantilePositions",
signature(object = "CAGEset"),
function (object, clusters, qLow = 0.1, qUp = 0.9, useMulticore = FALSE, nrCores = NULL){

	pt <- .Platform$OS.type
	if(useMulticore == TRUE){
		if(pt == "unix"){
			if("parallel" %in% rownames(installed.packages()) == FALSE){
				stop("Cannot use multicore because package 'parallel' is not installed!")
			}
		}else{
			useMulticore = FALSE
			warning("Multicore is not supported on non-Unix platforms! Setting useMulticore=FALSE")
		}
	}
	
	objName <- deparse(substitute(object))
	sample.labels = sampleLabels(object)
	
	message("\nGetting positions of quantiles within clusters...")
	
	ctss.clusters.q.low.list <- list()
	ctss.clusters.q.up.list <- list()
	
	if(clusters == "tagClusters"){	
		
		samples.cumsum.list <- object@CTSScumulativesTagClusters
		
		for(s in sample.labels) {
			
			message("\t-> ", s)
			
			clusters.cumsum.list <- samples.cumsum.list[[s]]
			ctss.clusters <- tagClusters(object, sample = s)
			ctss.clusters.q.low <- .get.quant.pos(cluster.cumsums = clusters.cumsum.list, coors = ctss.clusters, q = qLow, q.orientation = "low", use.multicore = useMulticore, nrCores = nrCores)
			ctss.clusters.q.up <- .get.quant.pos(cluster.cumsums = clusters.cumsum.list, coors = ctss.clusters, q = qUp, q.orientation = "up", use.multicore = useMulticore, nrCores = nrCores)

			ctss.clusters.q.low.list[[s]] <- ctss.clusters.q.low[,c(which(colnames(ctss.clusters.q.low) == "cluster"), grep("q_", colnames(ctss.clusters.q.low), fixed = T))]
			ctss.clusters.q.up.list[[s]] <- ctss.clusters.q.up[,c(which(colnames(ctss.clusters.q.up) == "cluster"), grep("q_", colnames(ctss.clusters.q.up), fixed = T))]
		
		}
		
		object@tagClustersQuantileLow <- ctss.clusters.q.low.list
		object@tagClustersQuantileUp <- ctss.clusters.q.up.list
		
	}else if (clusters == "consensusClusters"){
		
		samples.cumsum.list <- object@CTSScumulativesConsensusClusters
		ctss.clusters.orig <- merge(object@consensusClusters, object@consensusClustersTpmMatrix, by.x = 1, by.y = 0)
		
		for(s in sample.labels) {
			
			message("\t-> ", s)
			
			clusters.cumsum.list <- samples.cumsum.list[[s]]
			ctss.clusters <- ctss.clusters.orig[ctss.clusters.orig[,s]>0,]
			colnames(ctss.clusters)[which(colnames(ctss.clusters) == "consensus.cluster")] = "cluster"
			ctss.clusters.q.low <- .get.quant.pos(cluster.cumsums = clusters.cumsum.list, coors = ctss.clusters, q = qLow, q.orientation = "low", use.multicore = useMulticore, nrCores = nrCores)
			ctss.clusters.q.up <- .get.quant.pos(cluster.cumsums = clusters.cumsum.list, coors = ctss.clusters, q = qUp, q.orientation = "up", use.multicore = useMulticore, nrCores = nrCores)
			
			ctss.clusters.q.low.list[[s]] <- ctss.clusters.q.low[,c(which(colnames(ctss.clusters.q.low) == "cluster"), grep("q_", colnames(ctss.clusters.q.low), fixed = T))]
			ctss.clusters.q.up.list[[s]] <- ctss.clusters.q.up[,c(which(colnames(ctss.clusters.q.up) == "cluster"), grep("q_", colnames(ctss.clusters.q.up), fixed = T))]
			
		}
		
		object@consensusClustersQuantileLow <- ctss.clusters.q.low.list
		object@consensusClustersQuantileUp <- ctss.clusters.q.up.list
		
		
		
	}else{
		stop("'clusters' parameter must be one of the (\"tagClusters\", \"consensusClusters\")")
	}
	
	assign(objName, object, envir = parent.frame())
	invisible(1)
	

}
)
