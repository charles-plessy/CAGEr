#####
# Function that calculates reversed cumulative sums given a list of cumulative sums
# ARGUMENTS: cumsum.list - list of Rle vectors (IRanges package) with cumulative sums (first number in the vector needs to be a zero) (such as returned by 'get.cumsum' function)
# RETURNS: list of Rle vectors containing reversed cumulative sums for all elements in the input cumsum list (Rle vectors are shorter by 1 than original vectors because first zero (i.e. here last number) is omitted) 

.reverse.cumsum <- function(cumsum.list, use.multicore = F, nrCores = NULL) {
	
	if(use.multicore){
		library(multicore)
		if(is.null(nrCores)){
			nrCores <- multicore:::detectCores()
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
		library(multicore)
		if(is.null(nrCores)){
			nrCores <- multicore:::detectCores()
		}				
		scores <- unlist(mclapply(cum.sum.stages, function(x) {less.tpm <- which(x[nrow(x),] == min(x[nrow(x),]))[1]; if(max(x[,less.tpm])>0) {max(x[,less.tpm] - x[,(3-less.tpm)])/max(x[,less.tpm])}else{return(as.numeric(NA))}}, mc.cores = nrCores))
	}else{
		scores <- unlist(lapply(cum.sum.stages, function(x) {less.tpm <- which(x[nrow(x),] == min(x[nrow(x),]))[1]; if(max(x[,less.tpm])>0) {max(x[,less.tpm] - x[,(3-less.tpm)])/max(x[,less.tpm])}else{return(as.numeric(NA))}}))	
	}

}	
	



.getCAGEtagCountCoverage <- function(ctss, coors) {
	
	cov = rep(0, max(coors$end) + 1)
	cov[ctss$pos] = ctss$tagcount
	cov <- Rle(cov)
	cluster.tags <- Views(cov, start = coors$start, end = coors$end)
	cluster.tags <- viewSums(cluster.tags)
	return(as.integer(cluster.tags))
		
}


.getCAGEtagCountChr <- function(ctss.clusters, ctss.df, chrom, str) {
	
	clusters.chr.strand.coor <- subset(ctss.clusters, chr == chrom & strand == str)
	ctss.chr.strand <- subset(ctss.df, chr == chrom & strand == str)
	if(nrow(clusters.chr.strand.coor)>0){
		strand.tagcount <- .getCAGEtagCountCoverage(ctss = ctss.chr.strand, coors = clusters.chr.strand.coor)
		return(strand.tagcount)
	}else{
		return()
	}
	
}


#####
# Function that calculates total tag count in CAGE clusters
# ARGUMENTS: ctss.df - data frame with one row per CTSS containing at least five columns, *cluster (cluster ID) *chr (chromosome) *pos (genomic position of CTSSs) *strand (genomic strand) *tagcount (raw CAGE tag count)
#            ctss.clusters - data frame with one row per cluster containing at least 6 columns, *cluster (cluster ID) *chr (chromosome) *start (start position of the cluster) *end (end position of the cluster) *strand (strand) *dominant_ctss (position of dominant peak)
# RETURNS: integer vector of total tag count per cluster 


.getTotalTagCount <- function(ctss.df, ctss.clusters, id.column, use.multicore = FALSE, nrCores = NULL) {
	
	if(use.multicore == TRUE) {
		library(multicore)
		if(is.null(nrCores)){
			nrCores <- multicore:::detectCores()
		}		
		
		clusters.tagcount <- mclapply(as.list(unique(ctss.clusters$chr)), function(x) {
									
									plus.tagcount <- .getCAGEtagCountChr(ctss.clusters = ctss.clusters, ctss.df = ctss.df, chrom = x, str = "+")
									minus.tagcount <- .getCAGEtagCountChr(ctss.clusters = ctss.clusters, ctss.df = ctss.df, chrom = x, str = "-")
									tagcount.chr <- append(plus.tagcount, minus.tagcount)
									return(tagcount.chr)
									
									}, mc.cores = nrCores
									)
		clusters.tagcount <- unlist(clusters.tagcount)
		
	}else{
		
		clusters.tagcount <- list()
		
		for(chrom in unique(ctss.clusters$chr)) {
			
			plus.tagcount <- .getCAGEtagCountChr(ctss.clusters = ctss.clusters, ctss.df = ctss.df, chrom = chrom, str = "+")
			minus.tagcount <- .getCAGEtagCountChr(ctss.clusters = ctss.clusters, ctss.df = ctss.df, chrom = chrom, str = "-")
			clusters.tagcount <- append(clusters.tagcount, append(plus.tagcount, minus.tagcount))
			
		}
	}
	
	n <- unlist(sapply(unique(ctss.clusters$chr), function(x) {a = c(subset(ctss.clusters, chr == x & strand == '+')[,id.column], subset(ctss.clusters, chr == x & strand == '-')[,id.column]); return(a)}))
	names(clusters.tagcount) <- n
	return(clusters.tagcount)
	
}





.ksStat <- function(cumsum.matrix){
	
	z <- cumsum.matrix[,1]/max(cumsum.matrix[,1]) - cumsum.matrix[,2]/max(cumsum.matrix[,2])
	ks.stat <- max(abs(z))
	return(ks.stat)
	
}


.pkstwo <- function(x, tol = 1e-06) {

	if (is.numeric(x)) 
		x <- as.double(x)
	else stop("argument 'x' must be numeric")
	p <- rep(0, length(x))
	p[is.na(x)] <- NA
	IND <- which(!is.na(x) & (x > 0))
	if (length(IND)){
		if(exists("C_pkstwo", envir=asNamespace("stats"))){
			p[IND] <- .C(stats:::C_pkstwo, length(x[IND]), p = x[IND], as.double(tol))$p
		}else if(exists("C_pKS2", envir=asNamespace("stats"))){
			p[IND] <- .Call(stats:::C_pKS2, p = x[IND], as.double(tol))
		}
	
	}
	return(p)

}


.ksPvalue <- function(d, n){

	1 - .pkstwo(sqrt(n) * d)

}



