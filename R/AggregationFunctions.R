#' @include CAGEr.R

###################################################################
# Functions for aggregating tag clusters (TCs) across all samples
#


.make.consensus.clusters <- function(TC.list, plus.minus = 0, tpm.th = 0) {
	
	TC.df <- do.call(rbind, TC.list)
	TC.df <- subset(TC.df, tpm >= tpm.th)
	TC.gr <- GRanges(seqnames = TC.df$chr, ranges = IRanges(start = pmax(1, TC.df$start - plus.minus), end = TC.df$end + plus.minus), strand = TC.df$strand)
	consensus.clusters.gr <- reduce(TC.gr)
	
	consensus.clusters = data.frame()
	for(i in 1:length(TC.list)) {
		
		clusters.q = subset(TC.list[[i]], tpm >= tpm.th)
		col.n = colnames(clusters.q)
		if(nrow(clusters.q) > 0) {
			
			clusters.q.gr = GRanges(seqnames = clusters.q$chr, ranges = IRanges(start = pmax(1, clusters.q$start - plus.minus), end = clusters.q$end + plus.minus), strand = clusters.q$strand, values = clusters.q)
			o = findOverlaps(consensus.clusters.gr, clusters.q.gr)
			consensus.clusters = rbind(consensus.clusters, data.frame(consensus.cluster = queryHits(o), as.data.frame(values(clusters.q.gr[as.vector(subjectHits(o))])), sample = names(TC.list)[i]))
			
		}
		
	}
	
	colnames(consensus.clusters) = c('consensus.cluster', col.n, 'sample')
	return(consensus.clusters[order(consensus.clusters$consensus.cluster),])
	
}

#' @name consensusClusterConvertors
#' 
#' @title Private functions to convert CC formats
#' 
#' @details Interconvert consensus clusters (CC) formats used in classes CAGEset
#' (\code{data.frame}) and CAGEexp (\code{GRanges}).
#' 
#' @examples 
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' df <- consensusClusters(exampleCAGEset)
#' head(df)
#' gr <- CCdataframe2granges(df)
#' gr
#' # No round-trip because start and end were not integer in df.
#' # identical(df, CCgranges2dataframe(gr))
NULL

#' @name CCgranges2dataframe
#' 
#' @rdname consensusClusterConvertors
#' 
#' @param gr Consensus clusters in \code{GRanges} format.
#' 
#' @importFrom S4Vectors mcols
#' 
#' @export

CCgranges2dataframe <- function(gr) {
  gr$score <- NULL
  df <- data.frame( consensus.cluster  = gr$consensus.cluster)
  gr$consensus.cluster <- NULL
  if(!is.null(df$cluster)) {
    df <- cbind(df, gr$cluster)
    gr$cluster <- NULL
  }
  df <- cbind(df , data.frame( chr     = decode(seqnames(gr))
                             , start   = start(gr)
                             , end     = end(gr)
                             , strand  = decode(droplevels(strand(gr)))))
  df <- cbind(df, mcols(gr))
  as.data.frame(df)
}

#' @name CCdataframe2granges
#' @rdname consensusClusterConvertors
#' 
#' @param df Consensus clusters in \code{data.frame} format.
#' 
#' @family df2granges converters
#' 
#' @export

CCdataframe2granges <- function(df) {
	gr <- GRanges( seqnames           = df$chr
	             , ranges             = IRanges(df$start, df$end)
 	             , score              = df$tpm
               , strand             = df$strand)
	mcols(gr) <- cbind( mcols(gr)
	                  , df[,setdiff(colnames(df), c("chr", "start", "end", "strand")), drop = FALSE])
	names(gr) <- rownames(df)
	gr
}
