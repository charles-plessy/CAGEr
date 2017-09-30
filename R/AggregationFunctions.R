#' @include CAGEr.R

###################################################################
# Functions for aggregating tag clusters (TCs) across all samples
#

#' @name .make.consensus.clusters
#' @noRd
#' @importFrom IRanges reduce
#' @importFrom S4Vectors endoapply

setGeneric( ".make.consensus.clusters"
          , function(TC.list, plus.minus = 0, tpm.th = 0)
              standardGeneric(".make.consensus.clusters"))

setMethod(".make.consensus.clusters", "GRangesList", function(TC.list, plus.minus, tpm.th) {
  gr.list <- endoapply(TC.list, function (gr) gr <- gr[score(gr) >= tpm.th])
  
  clusters.gr <- unlist(gr.list)
  mcols(clusters.gr) <- NULL
	names(clusters.gr) <- NULL
  suppressWarnings(start(clusters.gr) <- start(clusters.gr) - plus.minus) # Suppress warnings because we trim later
	suppressWarnings(end(clusters.gr)   <- end(clusters.gr)   + plus.minus)
	clusters.gr <- reduce(trim(clusters.gr))
	
	gr.list <- endoapply( gr.list
                      , function (gr) {
    o = findOverlaps(clusters.gr, gr)
    gr$consensus.cluster <- queryHits(o)
    gr$sample <- "tbd" # Can not retreive 
    gr
  })
	
	for (i in seq_along(gr.list)) gr.list[[i]]$sample <- names(gr.list)[[i]]
	
	unname(unlist(gr.list))
})

#' @name consensusClusterConvertors
#' 
#' @title Private functions to convert CC formats
#' 
#' @description  Interconvert consensus clusters (CC) formats used in classes CAGEset
#' (\code{data.frame}) and CAGEexp (\code{GRanges}).
#' 
#' @examples 
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
  if (is.null(gr$tpm)) gr$tpm <- gr$score
  gr$score <- NULL
  consensus.cluster <-
    suppressWarnings(as.integer(gr$consensus.cluster))  # Make sure it does not sort lexically!
  if (any(is.na(consensus.cluster))) {                  # Revert if the IDs were not numbers.
    consensus.cluster <- gr$consensus.cluster
  }
  df <- data.frame(consensus.cluster = consensus.cluster) 
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
