setGeneric(
name="cumulativeCTSSdistribution",
def=function(object, clusters, useMulticore = FALSE, nrCores = NULL){
	standardGeneric("cumulativeCTSSdistribution")
}
)

setMethod("cumulativeCTSSdistribution",
signature(object = "CAGEset"),
function (object, clusters, useMulticore = FALSE, nrCores = NULL){

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

	message("\nCalculating cumulative sum of CAGE signal along clusters...")
	
	idx <- object@filteredCTSSidx
	ctss.df <- cbind(CTSScoordinates(object)[idx,], object@normalizedTpmMatrix[idx,,drop=F])

	samples.cumsum.list <- list()

	if(clusters == "tagClusters"){	
		
		for(s in sample.labels) {
		
			message("\t-> ", s)
			d <- ctss.df[,c("chr", "pos", "strand", s)]
			colnames(d) <- c("chr", "pos", "strand", "tpm")
            #d <- subset(d, tpm>0)
			ctss.clusters <- tagClusters(object, sample = s)
			clusters.cumsum.list <- .getCumsum(ctss.df = d, ctss.clusters = ctss.clusters, id.column = "cluster", use.multicore = useMulticore, nrCores = nrCores)
			samples.cumsum.list[[s]] <- clusters.cumsum.list
			invisible(gc())
						
		}
		
		object@CTSScumulativesTagClusters <- samples.cumsum.list
		
	}else if (clusters == "consensusClusters"){

		ctss.clusters.orig <- merge(object@consensusClusters, object@consensusClustersTpmMatrix, by.x = 1, by.y = 0)

		for(s in sample.labels) {
			
			message("\t-> ", s)
			d <- ctss.df[,c("chr", "pos", "strand", s)]
			colnames(d) <- c("chr", "pos", "strand", "tpm")
			d <- subset(d, tpm>0)
			ctss.clusters <- ctss.clusters.orig[ctss.clusters.orig[,s]>0,]
			clusters.cumsum.list <- .getCumsum(ctss.df = d, ctss.clusters = ctss.clusters, id.column = "consensus.cluster", use.multicore = useMulticore, nrCores = nrCores)
			samples.cumsum.list[[s]] <- clusters.cumsum.list
			invisible(gc())
			
		}
		
		object@CTSScumulativesConsensusClusters <- samples.cumsum.list
		
	}else{
		stop("'clusters' parameter must be one of the (\"tagClusters\", \"consensusClusters\")")
	}
	
	assign(objName, object, envir = parent.frame())
	invisible(1)
	
}
)
