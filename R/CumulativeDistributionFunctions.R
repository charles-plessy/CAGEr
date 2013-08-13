.getCAGEsignalCoverage <- function(ctss, coors) {
	
	cov = rep(0, max(coors$end) + 1)
	cov[ctss$pos+1] = ctss$tpm
	cov = cumsum(cov)
	return(cov)
	
}

.getCumsumChr <- function(cov.cumsum, coors) {
	
	cov.cumsum <- Rle(cov.cumsum)
	if(nrow(coors)>1) {
		cluster.cumsums <- Views(cov.cumsum, start = coors$start+1, end = coors$end+1)
		cluster.cumsums <- viewApply(cluster.cumsums, FUN = function(x) {x - x[1]})
	}else{
		cluster.cumsums <- list(cov.cumsum[(coors$start+1):(coors$end+1)] - cov.cumsum[coors$start+1])
	}
	return(cluster.cumsums)
	
}

.getCumsumChr2 <- function(ctss.clusters, ctss.df, chrom, str) {
	
	clusters.chr.strand.coor <- subset(ctss.clusters, chr == chrom & strand == str)
	ctss.chr.strand <- subset(ctss.df, chr == chrom & strand == str)
	if(nrow(clusters.chr.strand.coor)>0){
		cov.strand <- .getCAGEsignalCoverage(ctss = ctss.chr.strand, coors = clusters.chr.strand.coor)
		strand.cumsum <- .getCumsumChr(cov.cumsum = cov.strand, coors = clusters.chr.strand.coor)	
		return(strand.cumsum)
	}else{
		return()
	}
	
}


#####
# Function that calculates cumulative sums of tpm along the clusters
# ARGUMENTS: ctss.df - data frame with one row per CTSS containing at least five columns, *cluster (cluster ID) *chr (chromosome) *pos (genomic position of CTSSs) *strand (genomic strand) *tpm (CAGE tag count or number per million)
#            ctss.clusters - data frame with one row per cluster containing at least 6 columns, *cluster (cluster ID) *chr (chromosome) *start (start position of the cluster) *end (end position of the cluster) *strand (strand) *dominant_ctss (position of dominant peak)
# RETURNS: list of Rle vectors (IRanges package) containing cumulative sum for each cluster (length of list is equal to number of clusters and names of the list components corespond to the name of the corresponding cluster) v


.getCumsum <- function(ctss.df, ctss.clusters, id.column, use.multicore = FALSE, nrCores = NULL) {
		
	if(use.multicore == TRUE) {
		library(parallel)
		if(is.null(nrCores)){
			nrCores <- detectCores()
		}		
		
		clusters.cumsum <- mclapply(as.list(unique(ctss.clusters$chr)), function(x) {
										
										plus.cumsum <- .getCumsumChr2(ctss.clusters = ctss.clusters, ctss.df = ctss.df, chrom = x, str = "+")
										minus.cumsum <- .getCumsumChr2(ctss.clusters = ctss.clusters, ctss.df = ctss.df, chrom = x, str = "-")
										cumsum.chr <- append(plus.cumsum, minus.cumsum)
										return(cumsum.chr)
									 
									 }, mc.cores = nrCores
									)
		clusters.cumsum <- unlist(clusters.cumsum)
		
	}else{
		
		clusters.cumsum <- list()
	
		for(chrom in unique(ctss.clusters$chr)) {
				
			plus.cumsum <- .getCumsumChr2(ctss.clusters = ctss.clusters, ctss.df = ctss.df, chrom = chrom, str = "+")
			minus.cumsum <- .getCumsumChr2(ctss.clusters = ctss.clusters, ctss.df = ctss.df, chrom = chrom, str = "-")
			clusters.cumsum <- append(clusters.cumsum, append(plus.cumsum, minus.cumsum))
			
		}
	}
	
	n <- unlist(sapply(unique(ctss.clusters$chr), function(x) {a = c(subset(ctss.clusters, chr == x & strand == '+')[,id.column], subset(ctss.clusters, chr == x & strand == '-')[,id.column]); return(a)}))
	names(clusters.cumsum) <- n
	return(clusters.cumsum)
	
}

