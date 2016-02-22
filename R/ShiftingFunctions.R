#####
# Function that calculates reversed cumulative sums given a list of cumulative sums
# ARGUMENTS: cumsum.list - list of Rle vectors (IRanges package) with cumulative sums (first number in the vector needs to be a zero) (such as returned by 'get.cumsum' function)
# RETURNS: list of Rle vectors containing reversed cumulative sums for all elements in the input cumsum list (Rle vectors are shorter by 1 than original vectors because first zero (i.e. here last number) is omitted) 

.reverse.cumsum <- function(cumsum.list, use.multicore = F, nrCores = NULL) {
	
	if(use.multicore){
		library(parallel)
		if(is.null(nrCores)){
			nrCores <- detectCores()
		}		
		cumsum.reversed <- mclapply(cumsum.list, function(x) {cumsum(rev(x[2:length(x)] - x[1:(length(x)-1)]))}, mc.cores = nrCores)
	}else{
		cumsum.reversed <- lapply(cumsum.list, function(x) {cumsum(rev(x[2:length(x)] - x[1:(length(x)-1)]))})
	}
	return(cumsum.reversed)
	
}


.get.dominant.ctss <- function(v, isCumulative = FALSE){
	if(all(unique(v) == 0)){
		idx <- NA
	}else{
		if(isCumulative){
			raw <- v[2:length(v)] - v[1:(length(v)-1)]
		}else{
			raw <- v
		}
	idx <- which(raw == max(raw))
	if(length(idx > 1)){
		idx <- idx[ceiling(length(idx)/2)]
	}
	}
	return(as.integer(idx))
}


.score.promoter.shifting <- function(cum.sum.stages, use.multicore = F, nrCores = NULL) {	
	
	if(use.multicore){
		library(parallel)
		if(is.null(nrCores)){
			nrCores <- detectCores()
		}				
		scores <- unlist(mclapply(cum.sum.stages, function(x) {less.tpm <- which(x[nrow(x),] == min(x[nrow(x),]))[1]; if(max(x[,less.tpm])>0) {max(x[,less.tpm] - x[,(3-less.tpm)])/max(x[,less.tpm])}else{return(as.numeric(NA))}}, mc.cores = nrCores))
	}else{
		scores <- unlist(lapply(cum.sum.stages, function(x) {less.tpm <- which(x[nrow(x),] == min(x[nrow(x),]))[1]; if(max(x[,less.tpm])>0) {max(x[,less.tpm] - x[,(3-less.tpm)])/max(x[,less.tpm])}else{return(as.numeric(NA))}}))	
	}

}	
	


#####
# Function that calculates total tag count in CAGE clusters
# ARGUMENTS: ctss.df - data frame with one row per CTSS containing at least four columns, *chr (chromosome) *pos (genomic position of CTSSs) *strand (genomic strand) *tagcount (raw CAGE tag count)
#            ctss.clusters - data frame with one row per cluster containing at least 6 columns, *cluster (cluster ID) *chr (chromosome) *start (start position of the cluster) *end (end position of the cluster) *strand (strand) *dominant_ctss (position of dominant peak)
# RETURNS: integer vector of total tag count per cluster 



.getTotalTagCount <- function(ctss.df, ctss.clusters, id.column) {
	
	ctss.clusters.gr <- GRanges(seqnames = ctss.clusters$chr, ranges = IRanges(start = ctss.clusters$start + 1, end = ctss.clusters$end), strand = ctss.clusters$strand, consensus.cluster = ctss.clusters$consensus.cluster)
	ctss.gr <- GRanges(seqnames = ctss.df$chr, ranges = IRanges(start = ctss.df$pos, end = ctss.df$pos), strand = ctss.df$strand, tpm = ctss.df$tagcount)
	o <- findOverlaps(ctss.clusters.gr, ctss.gr)
	tpm.dt <- data.table(tpm = ctss.gr$tpm[subjectHits(o)], consensus.cluster = ctss.clusters.gr$consensus.cluster[queryHits(o)])
	tpm.df <- as.data.frame(tpm.dt[,sum(tpm),by=consensus.cluster])
	cc.zero <- subset(ctss.clusters, !(consensus.cluster %in% tpm.df$consensus.cluster))$consensus.cluster
	tpm.df <- rbind(tpm.df, data.frame(consensus.cluster = cc.zero, V1 = rep(0, length(cc.zero))))
	tpm.df <- tpm.df[order(tpm.df$consensus.cluster),]
	tag.count <- tpm.df$V1
	names(tag.count) <- tpm.df$consensus.cluster
	return(tag.count)
	
}



.ksStat <- function(cumsum.matrix){
	
	z <- cumsum.matrix[,1]/max(cumsum.matrix[,1]) - cumsum.matrix[,2]/max(cumsum.matrix[,2])
	ks.stat <- max(abs(z))
	return(ks.stat)
	
}



.pKS2 <- function(x, tol){
	
	k_max <- sqrt(2-log(tol))
	
	for(i in c(1:length(x))){
		
		if(x[i] < 1){
		
			z <- as.double(-1 * pi^2/(8*x[i]^2))
			w <- as.double(log(x[i]))
			k <- seq(1, k_max - 10^(-10), 2)
			s <- as.double(sum(exp(k^2 * z - w)))
			x[i] <- as.double(s*sqrt(2*pi))

		}else{
		
			z <- -2 * x[i]^2
			s <- -1
			k <- 1
			old <- 0
			new <- 1
			while(abs(old - new) > tol){
				old <- new
				new <- new + 2 * s * exp(z * k^2)
				s <- -1 * s
				k <- k + 1
			}
			x[i] <- new
			
		}
		
	}
	
	return(x)
	
}

	
.pkstwo <- function(x, tol = 1e-06) {

	if (is.numeric(x)){ 
		x <- as.double(x)
	}else{
		stop("argument 'x' must be numeric")
	}
	p <- rep(0, length(x))
	p[is.na(x)] <- NA
	IND <- which(!is.na(x) & (x > 0))
	if (length(IND) > 0){
		p[IND] <- .pKS2(x = x[IND], tol = as.double(tol))
	}
	return(p)

}


.ksPvalue <- function(d, n){

	1 - .pkstwo(sqrt(n) * d)

}



