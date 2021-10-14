#' @include CAGEr.R

###################################################################
# Functions for aggregating tag clusters (TCs) across all samples
#

#' @name .make.consensus.clusters
#' @noRd
#' @importFrom IRanges reduce
#' @importFrom S4Vectors endoapply

setGeneric( ".make.consensus.clusters"
          , function(TC.list, fix.at = "iqw", qLow = NULL, qUp = NULL, plus.minus = 0, tpm.th = 0)#, exclSigBel.th = TRUE)
              standardGeneric(".make.consensus.clusters"))

setMethod(".make.consensus.clusters", "GRangesList", 
          function(TC.list, fix.at, qLow, qUp, plus.minus, tpm.th
                   # , exclSigBel.th
                   ) {
  # Filter out TCs with too low score.
  
  # Shouldn't this condition be >= instead of > ? The help mentions (verbatim)
  # --
  # First, TCs with signal >= tpmThreshold from all CAGE datasets are selected, 
  # and their 5' and 3' boundaries are determined based on provided qLow and 
  # qUp parameter (or the start and end coordinates, if they are set to NULL). 
  # --
  # Now changed from > to >=
  # 
  gr.list <- endoapply(TC.list, function (gr) gr <- gr[score(gr) >= tpm.th])
  
  # Aggregate clusters by expanding and merging TCs from all samples.
  clusters.gr <- unlist(gr.list)
  # should I put actually this in a function, since I need to use it later again, shiting by quantile positions?!
  if (fix.at == "iqw"){
    end(clusters.gr)   <- as.integer(start(clusters.gr) +
                                       mcols(clusters.gr)[[paste0("q_", qUp)]] - 1)
    start(clusters.gr) <- as.integer(start(clusters.gr) +
                                       mcols(clusters.gr)[[paste0("q_", qLow)]] - 1)  
  suppressWarnings(start(clusters.gr) <- start(clusters.gr)  - plus.minus) # Suppress warnings
  suppressWarnings(end(clusters.gr)   <- end(clusters.gr) + plus.minus) # because we trim later
  
  } else if (fix.at == "start"){
    suppressWarnings(start(clusters.gr) <- start(clusters.gr) - plus.minus) # Suppress warnings
    suppressWarnings(end(clusters.gr)   <- end(clusters.gr)   + plus.minus) # because we trim later
    
  } else{
    error(
      "The consensus clusters must be fixed either on the start or quntile positions")
  }
  mcols(clusters.gr) <- NULL
  names(clusters.gr) <- NULL
  clusters.gr <- reduce(trim(clusters.gr))
  
  # Annotate TCs with ID of the aggregated cluster they intersect with.
  gr.list <- endoapply( gr.list, function (gr) {
    gro <- gr
    ## 
    ## (snikumbh) Not clear why this is needed here? Because this restriction of TCs 
    ## to quantiles as boundaries is already done above, right?
    if (fix.at == "iqw" ){
      end(gro)   <- as.integer(start(gro) + mcols(gro)[[paste0("q_", qUp)]] - 1)
      start(gro) <- as.integer(start(gro) + mcols(gro)[[paste0("q_", qLow)]] - 1)
    }
    o = findOverlaps(clusters.gr, gro)
    tmp <- gr[subjectHits(o)]
    tmp$consensus.cluster <- queryHits(o)
    #gr$consensus.cluster[subjectHits(o)] <- queryHits(o)
    tmp$sample <- "tbd" # Can not retrieve
    tmp
    # gr$consensus.cluster <- queryHits(o)
    # gr$sample <- "tbd" # Can not retreive
    # gr
  })

  
  ####
  # Add back the sample name.
  for (i in seq_along(gr.list)) gr.list[[i]]$sample <- names(gr.list)[[i]]

  unname(unlist(gr.list))
})

#' @title consensus Cluster convertors
#' 
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
