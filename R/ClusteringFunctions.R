#' @include CAGEexp.R CTSS.R
NULL

#' Summarise CTSSs included in clusters
#' 
#' @param ctss A [`CTSS`] object.
#' 
#' @param clusters A [`TagClusters`], [`ConsensusClusters`] or any other
#' object implementing the [`GRanges`] class.
#' 
#' @return The `clusters` object with a new `dominant_CTSS` metadata in `CTSS`
#' format reporting the genomic coordinate and expression score of most
#' highly expressed position in each cluster, plus a `nr_ctss` metadata reporting
#' the number of expressed CTSSs in each cluster.
#' 
#' @importFrom S4Vectors queryHits subjectHits runLength runValue
#' @export
#' 
#' @examples 
#' # See also benchmarks/dominant_ctss.md
#' (ctss <- CTSS( 'chr1', IRanges(start = 1:10, end = 1:10)
#'              , '+', score = c(1, 0, 0, 1, 2, 0, 2, 1, 0, 1)))
#' (clusters <- GRanges( 'chr1', IRanges(start = c(1,9)
#'                     , end = c(8,10)), '+')) |> as("TagClusters")
#' 
#' # The function assumes that all CTSSes have a score above zero
#' .ctss_summary_for_clusters(ctss[score(ctss)>0], clusters)
#' # If not the case, it will give incorrect nr_ctss and  fail to remove singletons
#' .ctss_summary_for_clusters(ctss, clusters)
#' 
#' # The function needs its output to be sorted and is not going to check it.
#' .ctss_summary_for_clusters(rev(ctss), clusters)
#' .ctss_summary_for_clusters(ctss, rev(clusters))
#' 
#' # Ties are resolved with 5' preference for both plus and minus strands.
#' # This may create a small bias.
#' ctss_minus <- ctss
#' strand(ctss_minus) <- '-'
#' clusters_minus <- clusters
#' strand(clusters_minus) <- '-'
#' .ctss_summary_for_clusters(ctss_minus, clusters_minus)

.ctss_summary_for_clusters <- function(ctss, clusters) {
  # Match the clusters and the CTSS
  o <- findOverlaps(clusters, ctss)

  # number of CTSS per cluster
  rl <- rle(queryHits(o))$length
  
  # Where each run starts in the CTSS ranges object
  cluster_start_idx <- cumsum(c(1, head(rl, -1)))
  
  # Scores sorted by run and position
  grouped_scores <- extractList(score(ctss), o)
  if (any(lengths(grouped_scores) == 0))
    stop("Some clusters do not overlap any CTSS.")
  
  # Find relative position of dominant CTSS in each run.
  # In case of ties, take the central one.  Note it might bias + and - strands
  # in a different way.  This is vectorised over all runs, because an sapply()
  # loop on the elements of a List is very slow.  The scores are decoded
  # because which() is not vectorised on RleList objects.
  decoded_scores <- extractList(decode(score(ctss)), o)
  max_pos <- which(decoded_scores == max(decoded_scores))
  local_max_idx <- unlist(max_pos)[ cumsum(c(0, head(lengths(max_pos), -1)))
                                  + ceiling(lengths(max_pos) / 2)]
  
  # Find absolute position of dominant CTSS in each run.
  global_max_ids <- cluster_start_idx + local_max_idx - 1
  
  # Record total expression of the cluster
  score(clusters) <- Rle(sum(grouped_scores))
  
  # Record dominant CTSS as CTSS object.
  clusters$dominant_ctss <- CTSS(granges(ctss)[subjectHits(o)][global_max_ids])
  
  # Record dominant CTSS score.
  score(clusters$dominant_ctss) <- score(ctss)[subjectHits(o)][global_max_ids]

  # Count the number of clusters   
  clusters$nr_ctss <- rl
  
  # Give numerical names to the clusters
  names(clusters) <- seq_along(clusters)
  
  # Finally, return the object
  clusters
}
