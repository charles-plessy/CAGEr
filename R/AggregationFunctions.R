###################################################################
# Functions for aggregating tag clusters (TCs) across all samples
#


.make.consensus.clusters <- function(TC.list, start.coor, end.coor, plus.minus = 0, tpm.th = 0) {
	
	TC.df <- do.call(rbind, TC.list)
	TC.df <- subset(TC.df, tpm >= tpm.th)
	TC.gr <- GRanges(seqnames = TC.df$chr, ranges = IRanges(start = pmax(1, TC.df[,start.coor] - plus.minus), end = TC.df[,end.coor] + plus.minus), strand = TC.df$strand)
	consensus.clusters.gr <- reduce(TC.gr)
	
	consensus.clusters = data.frame()
	for(i in 1:length(TC.list)) {
		
		clusters.q = subset(TC.list[[i]], tpm >= tpm.th)
		col.n = colnames(clusters.q)
		if(nrow(clusters.q) > 0) {
			
			clusters.q.gr = GRanges(seqnames = clusters.q$chr, ranges = IRanges(start = pmax(1, clusters.q[,start.coor] - plus.minus), end = clusters.q[,end.coor] + plus.minus), strand = clusters.q$strand, values = clusters.q)
			o = findOverlaps(consensus.clusters.gr, clusters.q.gr)
			consensus.clusters = rbind(consensus.clusters, data.frame(consensus.cluster = queryHits(o), as.data.frame(values(clusters.q.gr[as.vector(subjectHits(o))])), sample = names(TC.list)[i]))
			
		}
		
	}
	
	colnames(consensus.clusters) = c('consensus.cluster', col.n, 'sample')
	return(consensus.clusters[order(consensus.clusters$consensus.cluster),])
	
}


