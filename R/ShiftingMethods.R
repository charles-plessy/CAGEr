setGeneric(
name="scoreShift",
def=function(object, groupX, groupY, useMulticore = F, nrCores = NULL){
	standardGeneric("scoreShift")
}
)

setMethod("scoreShift",
signature(object = "CAGEset", groupX = "character", groupY = "character"),
function (object, groupX, groupY, useMulticore = F, nrCores = NULL){
	
	pt <- .Platform$OS.type
	if(useMulticore == TRUE){
		if(pt == "unix"){
			if("multicore" %in% rownames(installed.packages()) == FALSE){
				stop("Cannot use multicore because package 'multicore' is not installed!")
			}
		}else{
			useMulticore = FALSE
			warning("Multicore is not supported on non-Unix platforms! Setting useMulticore=FALSE")
		}
	}
	
	objName <- deparse(substitute(object))
	sample.labels <- sampleLabels(object)
	
	if(!(all(groupX %in% sample.labels) & all(groupY %in% sample.labels))){
		stop("'groupX' and 'groupY' must be vectors of sample labels! Check 'sampleLabels(", objName, ")' for labels of your CAGE datasets!")
	}
	
	a <- object@CTSScumulativesConsensusClusters
	b <- object@consensusClusters
	if(useMulticore){
		library(multicore)
		if(is.null(nrCores)){
			nrCores <- multicore:::detectCores()
		}		
		cumsum.list <- mclapply(a, function(x) {n <- names(x); y <- subset(b, !(consensus.cluster %in% as.integer(names(x)))); nulls <- lapply(as.list(c(1:nrow(y))), function(t) {Rle(rep(0, y[t, "end"] - y[t, "start"] + 1))}); x <- append(x, nulls); names(x) <- c(n, as.character(y$consensus.cluster)); return(x)}, mc.cores = nrCores)
	}else{
		cumsum.list <- lapply(a, function(x) {n <- names(x); y <- subset(b, !(consensus.cluster %in% as.integer(names(x)))); nulls <- lapply(as.list(c(1:nrow(y))), function(t) {Rle(rep(0, y[t, "end"] - y[t, "start"] + 1))}); x <- append(x, nulls); names(x) <- c(n, as.character(y$consensus.cluster)); return(x)})
	}
	cumsum.list.r <- lapply(cumsum.list, function(x) {.reverse.cumsum(x, use.multicore = useMulticore)})
	
	if(useMulticore){
		cumsum.matrices.list.f <- mclapply(as.list(names(cumsum.list[[1]])), function(x) {sapply(names(cumsum.list), function(y) {as.numeric(cumsum.list[[y]][[x]])})}, mc.cores = nrCores)
		cumsum.matrices.list.r <- mclapply(as.list(names(cumsum.list.r[[1]])), function(x) {m <- matrix(sapply(names(cumsum.list.r), function(y) {as.numeric(cumsum.list.r[[y]][[x]])}), ncol = length(names(cumsum.list.r))); colnames(m) <- names(cumsum.list.r); return(m)}, mc.cores = nrCores)
		cumsum.matrices.groups.f <- mclapply(cumsum.matrices.list.f, function(x) {cbind(groupX = rowSums(x[,groupX,drop=F]), groupY = rowSums(x[,groupY,drop=F]))}, mc.cores = nrCores)
		cumsum.matrices.groups.r <- mclapply(cumsum.matrices.list.r, function(x) {cbind(groupX = rowSums(x[,groupX,drop=F]), groupY = rowSums(x[,groupY,drop=F]))}, mc.cores = nrCores)
	}else{
		cumsum.matrices.list.f <- lapply(as.list(names(cumsum.list[[1]])), function(x) {sapply(names(cumsum.list), function(y) {as.numeric(cumsum.list[[y]][[x]])})})		
		cumsum.matrices.list.r <- lapply(as.list(names(cumsum.list.r[[1]])), function(x) {m <- matrix(sapply(names(cumsum.list.r), function(y) {as.numeric(cumsum.list.r[[y]][[x]])}), ncol = length(names(cumsum.list.r))); colnames(m) <- names(cumsum.list.r); return(m)})
		cumsum.matrices.groups.f <- lapply(cumsum.matrices.list.f, function(x) {cbind(groupX = rowSums(x[,groupX,drop=F]), groupY = rowSums(x[,groupY,drop=F]))})
		cumsum.matrices.groups.r <- lapply(cumsum.matrices.list.r, function(x) {cbind(groupX = rowSums(x[,groupX,drop=F]), groupY = rowSums(x[,groupY,drop=F]))})
	}
	
	names(cumsum.matrices.groups.f) <- names(cumsum.list[[1]])
	names(cumsum.matrices.groups.r) <- names(cumsum.list[[1]])
	
	if(useMulticore){
		dominant.ctss.pos <- mclapply(as.list(names(cumsum.matrices.groups.f)), function(x) {sapply(c("groupX", "groupY"), function(y) {.get.dominant.ctss(cumsum.matrices.groups.f[[x]][,y], isCumulative = T)})}, mc.cores = nrCores)
	}else{
		dominant.ctss.pos <- lapply(as.list(names(cumsum.matrices.groups.f)), function(x) {sapply(c("groupX", "groupY"), function(y) {.get.dominant.ctss(cumsum.matrices.groups.f[[x]][,y], isCumulative = T)})})	
	}
	dominant.ctss.pos <- data.frame(consensus.cluster = as.integer(names(cumsum.matrices.groups.f)), do.call(rbind, dominant.ctss.pos))
	dominant.ctss.pos <- dominant.ctss.pos[order(dominant.ctss.pos$consensus.cluster),]
	colnames(dominant.ctss.pos) <- c("consensus.cluster", "groupX.pos", "groupY.pos")

	clusters.info <- merge(b[,c(1:5)], dominant.ctss.pos, by.x = "consensus.cluster", by.y = "consensus.cluster")
	clusters.info$groupX.pos <- clusters.info$groupX.pos + clusters.info$start
	clusters.info$groupY.pos <- clusters.info$groupY.pos + clusters.info$start

	n <- names(cumsum.matrices.groups.f)
	cumsum.matrices.groups.f <- lapply(cumsum.matrices.groups.f, function(x) {x[-1,,drop=F]})
	
	scores.f <- .score.promoter.shifting(cumsum.matrices.groups.f, use.multicore = useMulticore, nrCores = nrCores)
	scores.r <- .score.promoter.shifting(cumsum.matrices.groups.r, use.multicore = useMulticore, nrCores = nrCores)
	
	scores <- pmax(scores.f, scores.r)
	names(scores) <- n

	groupX.tpm <- unlist(lapply(cumsum.matrices.groups.f, function(x) {max(x[,"groupX"])}))
	groupY.tpm <- unlist(lapply(cumsum.matrices.groups.f, function(x) {max(x[,"groupY"])}))
	scores.df <- data.frame(consensus.cluster = as.integer(names(scores)), shifting.score = scores, groupX.tpm = groupX.tpm, groupY.tpm = groupY.tpm)
	
	clusters.info <- merge(clusters.info, scores.df, by.x = "consensus.cluster", by.y = "consensus.cluster")
	clusters.info <- clusters.info[,c("consensus.cluster", "shifting.score", "groupX.pos", "groupY.pos", "groupX.tpm", "groupY.tpm")]
	
	object@consensusClustersShiftingScores <- clusters.info
	object@shiftingGroupX <- groupX
	object@shiftingGroupY <- groupY

	assign(objName, object, envir = parent.frame())
	invisible(1)

}
)


setGeneric(
name="getShiftingPromoters",
def=function(object, tpmThreshold, scoreThreshold){
	standardGeneric("getShiftingPromoters")
}
)

setMethod("getShiftingPromoters",
signature(object = "CAGEset"),
function (object, tpmThreshold, scoreThreshold){

	shifting.scores <- object@consensusClustersShiftingScores
	clusters <- object@consensusClusters
	
	sig.shifting <- subset(shifting.scores, (groupX.tpm >= tpmThreshold & groupY.tpm >= tpmThreshold) & shifting.score >= scoreThreshold)
	sig.shifting <- merge(clusters[,c(1:5)], sig.shifting, by.x = "consensus.cluster", by.y = "consensus.cluster", all.x = F, all.y = T)
	
	return(sig.shifting)
	
}
)

