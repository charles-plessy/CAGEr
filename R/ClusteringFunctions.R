#' @include CAGEexp.R CTSS.R

####################################################################################
# Implementation of simple distance-based clustering of data attached to sequences
#

#' @name distclu-functions
#' 
#' @title Private functions for distance clustering.
#' 
#' @param max.dist See [clusterCTSS()].
#' @param useMulticore,nrCores See clusterCTSS.
#' 
#' @description The flow of data is that a [CTSS] object of CTSSes is
#' progressively deconstructed, and data to form the clusters is progressively
#' integrated in a [data.table] object, which is finally converted to [GRanges]
#' at the end.  Doing the whole clustering with `GRanges` is more elegant, but
#' looping on a `GRangesList` was just too slow.  Maybe the operation on the
#' `data.table` is more efficient because it is vectorised.
#' 
#' @examples
#' # Get example data
#' library(IRanges)
#' library(GenomicRanges)
NULL

#' @name .cluster.ctss.strand
#' @rdname distclu-functions
#' 
#' @description `.cluster.ctss.strand` does the strandless distance clustering
#' of strandless CTSS positions from a single chromosome.  Input does not need
#' to be sorted, but _pay attention that the output is sorted_.
#' 
#' @param ctss.ipos.chr A IPos object.
#' 
#' @return `.cluster.ctss.strand` returns an [data.table] object containing
#' arbitrary cluster IDs (as integers) for each CTSS.
#'  
#' @importFrom data.table data.table
#' @importFrom IRanges IPos
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors runLength
#' @importFrom S4Vectors runValue
#' @importFrom S4Vectors subjectHits
#' 
#' @examples
#' 
#' #.cluster.ctss.strand
#' ctss.ipos.chr <- IPos(c(1,3,4,12,14,25,28))
#' CAGEr:::.cluster.ctss.strand(ctss.ipos.chr, 5)
#' 
#' ctss.chr <- CTSScoordinatesGR(exampleCAGEexp)
#' ctss.chr <- ctss.chr[strand(ctss.chr) == "+"]
#' ctss.ipos.chr <- ranges(ctss.chr)
#' # Same result if not sorted
#' identical(
#'   CAGEr:::.cluster.ctss.strand(ctss.ipos.chr, 20),
#'   CAGEr:::.cluster.ctss.strand(ctss.ipos.chr[sample(seq_along(ctss.ipos.chr))], 20)
#' )
#' # Returns an emtpy data.table object if given an empty IRanges object.
#' identical(CAGEr:::.cluster.ctss.strand(IPos(), 20), data.table::data.table())

setGeneric(".cluster.ctss.strand", function(ctss.ipos.chr, max.dist) standardGeneric(".cluster.ctss.strand"))

setMethod(".cluster.ctss.strand", "IPos", function(ctss.ipos.chr, max.dist) {
  if (identical(ctss.ipos.chr, IPos())) return(data.table())
  v <- Rle(rep(TRUE, max(start(ctss.ipos.chr)) + max.dist))
  v[start(ctss.ipos.chr) + max.dist] <- FALSE
  longRuns <- which(runLength(v) >= max.dist & runValue(v))
  c.starts <- cumsum(runLength(v))[longRuns] - max.dist
  c.ends <- c(cumsum(runLength(v))[(longRuns - 1)[-1]], length(v)) - max.dist
  clusters <- IRanges(start = c.starts, end = c.ends)
  o <- findOverlaps(clusters, ctss.ipos.chr)
  data.table(id = queryHits(o))
})


#' @name .cluster.ctss.chr
#' @rdname distclu-functions
#' 
#' @description \code{.cluster.ctss.chr} does the stranded distance clustering of CTSS on a
#' single chromosome, by dispatching both strands to \code{.cluster.ctss.strand} and merging
#' the results, taking care keep the cluster IDs unique.  Be careful that this function does
#' not look at the score.
#' 
#' @param ctss.chr A CTSS.chr object.
#' 
#' @return \code{.cluster.ctss.chr} returns a \code{\link{data.table}} object representing the
#' chromosome coordinates (\code{chr}, \code{pos}, \code{strand}) of each CTSS, with their
#' cluster ID (\code{id}).
#' 
#' @importFrom data.table data.table
#' @importFrom S4Vectors decode
#' 
#' @examples 
#' 
#' #.cluster.ctss.chr
#' ctss.chr <- as(CTSScoordinatesGR(exampleCAGEexp), "CTSS.chr")
#' CAGEr:::.cluster.ctss.chr(ctss.chr, 20)

setGeneric(".cluster.ctss.chr", function(ctss.chr, max.dist) standardGeneric(".cluster.ctss.chr"))

setMethod(".cluster.ctss.chr", "CTSS.chr", function(ctss.chr, max.dist) {
  clusters.p <- .cluster.ctss.strand(ranges(ctss.chr[strand(ctss.chr) == "+"]), max.dist)
  clusters.m <- .cluster.ctss.strand(ranges(ctss.chr[strand(ctss.chr) == "-"]), max.dist)
  maxIdP <- if (identical(clusters.p, data.table())) 0 else max(clusters.p)
  clusters.m <- clusters.m + maxIdP
  data.table( rbind(clusters.p, clusters.m)
            , chr    = as.character(seqnames(ctss.chr))
            , pos    = start(ctss.chr)
            , strand = decode(strand(ctss.chr)))
})


#' @name .ctss2clusters
#' @rdname distclu-functions
#' 
#' @description \code{.ctss2clusters} does the stranded distance clustering of CTSS.
#' 
#' @param ctss A CTSS object with a score column.
#' 
#' @return \code{.ctss2clusters} returns a \code{\link{data.table}} object representing the
#' cluster ID (\code{id}), chromosome coordinates (\code{chr}, \code{pos}, \code{strand}) and
#' the \code{score} of each CTSS.
#' 
#' @importFrom data.table data.table
#' @importFrom S4Vectors decode
#' 
#' @examples 
#' 
#' # .ctss2clusters
#' ctss <- CTSScoordinatesGR(exampleCAGEexp)
#' score(ctss) <- CTSSnormalizedTpmDF(exampleCAGEexp)[[1]]
#' seqnames(ctss)[rep(c(TRUE,FALSE), length(ctss) / 2)] <- "chr16"
#' ctss
#' clusters <- CAGEr:::.ctss2clusters(ctss, 20)
#' clusters

setGeneric(".ctss2clusters", function(ctss, max.dist = 20, useMulticore = FALSE, nrCores = NULL) standardGeneric(".ctss2clusters"))

setMethod(".ctss2clusters", "CTSS", function(ctss, max.dist, useMulticore, nrCores) {
  ctss <- sort(ctss)
  ctss <- ctss[score(ctss) != 0]
  ctss.list <- split(ctss, droplevels(seqnames(ctss)))
  ctss.list <- lapply(ctss.list, as, "CTSS.chr")
  ctss.list <- bplapply( ctss.list, .cluster.ctss.chr, max.dist = max.dist
                       , BPPARAM = CAGEr_Multicore(useMulticore, nrCores))
	max.clust <- sapply(ctss.list, function(x) {max(x$id)})
	max.clust <- cumsum(c(0, max.clust[-length(max.clust)]))
	
  for (x in seq_along(max.clust))
    ctss.list[[x]]$id <- ctss.list[[x]]$id + max.clust[x]
	
	dt <- do.call(rbind, ctss.list)
	dt$tpm <- decode(score(ctss))
	dt
})


#' @name .summarize.clusters
#' @rdname distclu-functions
#' @description \code{.summarize.clusters} calculates the number of CTSS, and
#' the position and score of a main peak, for each cluster.
#' 
#' @param ctss.clustered A \code{\link{data.table}} object representing the cluster ID
#'        (\code{id}), chromosome coordinates (\code{chr}, \code{pos}, \code{strand}) and
#'        the \code{score} of each CTSS.
#' @param removeSingletons Remove \dQuote{singleton} clusters that span only a single
#'        nucleotide ? (default = FALSE).
#' @param keepSingletonsAbove Even if \code{removeSingletons = TRUE}, keep singletons when
#'        their score is aboove threshold (default = \code{Inf}).
#' 
#' @return \code{.summarize.clusters} returns GRanges describing the clusters.
#' 
#' @importFrom data.table setnames
#' 
#' @examples 
#' 
#' # .summarize.clusters
#' CAGEr:::.summarize.clusters(clusters)
#' CAGEr:::.summarize.clusters(clusters, removeSingletons = TRUE)
#' CAGEr:::.summarize.clusters(clusters, removeSingletons = TRUE, keepSingletonsAbove = 5)

setGeneric(".summarize.clusters", function(ctss.clustered, removeSingletons = FALSE, keepSingletonsAbove = Inf) standardGeneric(".summarize.clusters"))

setMethod(".summarize.clusters", "data.table", function(ctss.clustered, removeSingletons, keepSingletonsAbove) {
	
  chr <- pos <- tpm <- cluster <- id <- NULL  # To keep R CMD check happy.
  clusters <- ctss.clustered[ , list( chr[1]
                                    , min(pos)
                                    , max(pos)
                                    , strand[1]
                                    , length(pos)
                                    , pos[find.dominant.idx(tpm)]
                                    , sum(tpm)
                                    , max(tpm))
                              , by = id]
  setnames(clusters, c( "cluster"
                      , "chr", "start", "end", "strand"
                      , "nr_ctss", "dominant_ctss", "tpm", "tpm.dominant_ctss"))
  
	if(removeSingletons)
		clusters <- subset(clusters, clusters$nr_ctss > 1 | clusters$tpm >= keepSingletonsAbove)
  
  gr <- GRanges( seqnames = Rle(factor(clusters$chr))
               , ranges   = IRanges(clusters$start, clusters$end)
               , strand   = clusters$strand
               , score    = Rle(clusters$tpm)
               , nr_ctss  = clusters$nr_ctss
               , dominant_ctss = clusters$dominant_ctss
               , tpm.dominant_ctss = Rle(clusters$tpm.dominant_ctss)
  )
  gr <- sort(gr)
  names(gr) <- seq_along(gr)
  gr
})


find.dominant.idx <- function (x) {
  w <- which(x == max(x))
  w[ceiling(length(w)/2)]
}


#' @name .distclu
#' @rdname distclu-functions
#' @description  \code{.distclu} receives the data from the main \code{clusterCTSS} and
#' dispatches each for (possibly parallel) processing.
#' 
#' @param se A \code{\link{SummarizedExperiment}} object representing the CTSSes and
#'        their expression in each sample.
#'         
#' @return \code{.distclu} returns GRanges describing the clusters.
#' 
#' @importFrom GenomicRanges GRangesList
#' @importFrom SummarizedExperiment rowRanges
#' 
#' @examples 
#' 
#' # .distclu
#' CAGEr:::.distclu(CTSStagCountSE(exampleCAGEexp))
#' \dontrun{
#' CAGEr:::.distclu(CTSStagCountSE(exampleCAGEexp), useMulticore = TRUE)
#' }

setGeneric(".distclu", function(se, max.dist = 20, removeSingletons = FALSE, keepSingletonsAbove = Inf, useMulticore = FALSE, nrCores = NULL) standardGeneric(".distclu"))

setMethod(".distclu", "SummarizedExperiment", function(se, max.dist, removeSingletons, keepSingletonsAbove, useMulticore, nrCores) {
	
  ctss.cluster.list <- list()
  for(s in colnames(se)) {
    message("\t-> ", s)
    d <- as(rowRanges(se), "CTSS")
    score(d) <- assays(se)[["normalizedTpmMatrix"]][[s]]
    d <- subset(d, score(d) > 0)
    clusters <- .ctss2clusters(ctss = d, max.dist = max.dist, useMulticore = useMulticore, nrCores = nrCores)
    ctss.cluster.list[[s]] <- clusters
  }
  ctss.cluster.list <- bplapply( ctss.cluster.list, .summarize.clusters
                               , removeSingletons = removeSingletons
                               , keepSingletonsAbove = keepSingletonsAbove
                               , BPPARAM = CAGEr_Multicore(useMulticore, nrCores))
  GRangesList(ctss.cluster.list)
})

#' @noRd
#' @examples 
#' # See also benchmarks/dominant_ctss.md
#' (ctss <- GRanges('chr1', IRanges(start = 1:10, end = 1:10), '+', score = c(1, 0, 0, 1, 2, 0, 2, 1, 0, 1)))
#' (clusters <- GRanges('chr1', IRanges(start = c(1,9), end = c(8,10)), '+'))
#' 
#' # The function assumes that all CTSSes have a score above zero
#' .ctss_summary_for_clusters(ctss[score(ctss)>0], clusters, removeSingletons = TRUE)
#' # If not the case, it will give incorrect nr_ctss and  fail to remove singletons
#' .ctss_summary_for_clusters(ctss, clusters, removeSingletons = TRUE)
#' 
#' # The function needs its output to be sorted and is not going to check it.
#' .ctss_summary_for_clusters(rev(ctss), clusters)
#' .ctss_summary_for_clusters(ctss, rev(clusters))
#' 
#' # Ties are resolved with 5' preference for both plus and minus strands.
#' # This may create a small bias.
#' .ctss_summary_for_clusters(ctss |> plyranges::mutate(strand = '-'), clusters |> plyranges::mutate(strand = '-'))

.ctss_summary_for_clusters <- function(ctss, clusters, removeSingletons = FALSE, keepSingletonsAbove = Inf) {
  # Match the clusters and the CTSS
  o <- findOverlaps(clusters, ctss)

  # number of CTSS per cluster
  rl <- rle(queryHits(o))$length
  
  # Where each run starts in the CTSS ranges object
  cluster_start_idx <- cumsum(c(1, head(rl, -1)))
  
  # Scores sorted by run and position
  grouped_scores <- extractList(score(ctss), o)
  
  # Find relative position of dominant CTSS in each run.
  # In case of ties, take the central one.  Note it might bias + and - strands
  # in a different way.
  local_max_idx <- sapply(grouped_scores, \(x) {
    w <- which(x == max(x))
    w[ceiling(length(w)/2)]
  })
  
  # Find absolute position of dominant CTSS in each run.
  global_max_ids <- cluster_start_idx + local_max_idx - 1
  
  # Record dominant CTSS as GRanges object.
  clusters$dominant_ctss <- granges(ctss)[subjectHits(o)][global_max_ids]
  
  # Record dominant CTSS score.  Mabye we should use its GRanges's score instead.
  clusters$tpm.dominant_ctss <-   score(ctss)[subjectHits(o)][global_max_ids]
  
  # Record total expression of the cluster
  score(clusters) <- Rle(sum(grouped_scores))

  # Count the number of clusters   
  clusters$nr_ctss <- rl
  
  # Remove clusters that match only one CTSS unless their expression is high enough
  if(removeSingletons)
    clusters <- subset(clusters, clusters$nr_ctss > 1 | score(clusters) >= keepSingletonsAbove)
  
  # Give numerical names to the clusters
  names(clusters) <- seq_along(clusters)
  
  # Finally, return the object
  clusters
}

setGeneric(".distclu2",
  function(
    object,
    max.dist = 20,
    removeSingletons = FALSE,
    keepSingletonsAbove = Inf,
    useMulticore = FALSE, nrCores = NULL
    ) standardGeneric(".distclu2"))

setMethod(".distclu2", "SummarizedExperiment",
          function(object, max.dist, removeSingletons, keepSingletonsAbove, useMulticore, nrCores) {
  ctss.cluster.list <- GRangesList()
  for(s in colnames(object)) {
    message("\t-> ", s)
    d <- as(rowRanges(object), "CTSS")
    score(d) <- assays(object)[["normalizedTpmMatrix"]][[s]]
    d <- subset(d, score(d) > 0)
    ctss.cluster.list[[s]] <-
      .distclu2(d, max.dist = max.dist,
                removeSingletons = removeSingletons,
                keepSingletonsAbove = keepSingletonsAbove)
  }
  ctss.cluster.list
})

.distclu2_CTSS <- function(object, max.dist, removeSingletons, keepSingletonsAbove, useMulticore, nrCores) {
  clusters <- reduce(GRanges(object), min = max.dist)
  clusters <- .ctss_summary_for_clusters(object, clusters,
                                         removeSingletons    = removeSingletons,
                                         keepSingletonsAbove = keepSingletonsAbove)
  names(clusters) <- seq_along(clusters)
  clusters
}

setMethod(".distclu2", "CTSS", .distclu2_CTSS)
