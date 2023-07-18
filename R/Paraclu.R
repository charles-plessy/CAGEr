#' Parametric clustering
#' 
#' Implementation of Paraclu - parametric clustering of data attached to
#' sequences (<http://www.cbrc.jp/paraclu/>).
#' 
#' @param object A [`CTSS`], or a [`S4Vectors::Pairs`] object with positions
#'        _first_ and scores _second_.
#' 
#' @param minStability Minimal stability of the cluster, where stability is
#'        defined as ratio between maximal and minimal density value for which
#'        this cluster is maximal scoring.  For definition of stability refer to
#'        Frith _et al._, Genome Research, 2007.  Clusters with stability
#'        `< minStability` will be discarded.
#' 
#' @param maxLength Maximal length of cluster in base-pairs.  Clusters with length
#'        `> maxLength` will be discarded.
#' 
#' @param removeSingletons Logical indicating if tag clusters containing only
#'        one CTSS be removed.
#' 
#' @param keepSingletonsAbove Controls which singleton tag clusters will be
#'        removed.  When `removeSingletons = TRUE`, only singletons with signal
#'        `< keepSingletonsAbove` will be removed.  Useful to prevent removing
#'        highly supported singleton tag clusters.  Default value `Inf` results
#'        in removing all singleton TCs when `removeSingletons = TRUE`.
#' 
#' @param reduceToNonoverlapping Logical, should smaller clusters contained
#'        within bigger cluster be removed to make a final set of tag clusters
#'        non-overlapping.
#' 
#' @param useMulticore Logical, should multicore be used.  `useMulticore = TRUE`
#'        has no effect on non-Unix-like platforms.
#' 
#' @param nrCores Number of cores to use when `useMulticore = TRUE`.  Default
#'        value `NULL` uses all detected cores.
#' 
#' @references 
#' MC Frith, E Valen, A Krogh, Y Hayashizaki, P Carninci, A Sandelin.  _A code
#' for transcription initiation in mammalian genomes._  Genome Research 2008
#' 18(1):1-12)
#' 
#' @returns Running Paraclu on a `Pairs` object containing positions and scores
#' returns an `IRanges` object containing the start and end positions of the
#' clusters, as well as the minimum and maximum density in `min_d` and `max_d`
#' metadata columns.
#' 
#' Running Paraclu on a `CTSS` object dispatches the computation on each strand
#' of each sequence level of the object, collects the `IRanges` and assemble
#' them back in a [`TagClusters`] object after filtering them by size and by
#' expression following the `minStability`, `maxLength`, `removeSingletons`,
#' `keepSingletonsAbove` and `reduceToNonoverlapping` parameters.
#' 
#' Running Paraclu on a [`RangedSummarizedExperiment`] object will loop on each
#' sample, and return the results as a [`GRangesList`] of `TagClusters`.
#' 
#' @family CAGEr clustering methods
#' 
#' @importFrom utils tail
#' 
#' @examples 
#' (ctss <- CTSSnormalizedTpmGR(exampleCAGEexp,1))
#' (pair <- Pairs(pos(ctss), score(ctss)))
#' CAGEr:::.paraclu_params(first(pair), second(pair))
#' CAGEr:::.paraclu(first(pair)[1:10], second(pair)[1:10])
#' paraclu(pair[1:10])
#' paraclu(ctss[1:10])
#' paraclu(CTSStagCountSE(exampleCAGEexp)[1:25,])
#'
#' @export

setGeneric("paraclu",
  function( object
          , minStability = 1, maxLength = 500
          , removeSingletons = FALSE, keepSingletonsAbove = Inf
          , reduceToNonoverlapping = TRUE
          , useMulticore = FALSE, nrCores = NULL)
   standardGeneric("paraclu"))

.paraclu_params <- function(pos, score) {
  sit   <- length(pos)
  tot   <- sum(score)
  if(sit == 1) return(list(br = NA, min_density = Inf, tot = tot, sit = sit))
  densities_forward <- cumsum(    score) [-sit] / (pos[2:sit] - pos[1]        )
  densities_reverse <- cumsum(rev(score))[-sit] / (pos[sit]   - pos[(sit-1):1])
  min_densities = c(min(densities_forward), min(densities_reverse))
  breaks <- c(     1 +   which(densities_forward == min_densities[1])[1]
                   , sit + 1 - which(densities_reverse == min_densities[2])[1])
  min_density <- min(min_densities)
  br <- breaks[tail(which(min_densities == min_density),1)]
  list(br = br, min_density = min_density, tot = tot, sit = sit)
}

.paraclu <- function( pos, score
                    , min_density = -Inf
                    , clusters = data.frame()) {
  params <- .paraclu_params(pos, score)
  
  if (!is.na(params$br)) {
    new_min     <- max(min_density, params$min_density)
    left  <-    1      : (params$br-1)
    right <- params$br :  length(pos)
    clusters <- rbind(.paraclu(pos[left],  score[left],  new_min, clusters), 
                      .paraclu(pos[right], score[right], new_min, clusters))
  }
  
  rbind( clusters
         , data.frame( start = min(pos)
                       ,   end = max(pos)
                       , min_d = min_density
                       , max_d = params$min_density))
}

#' @rdname paraclu

setMethod("paraclu", "Pairs",
  function( object
          , minStability = 1, maxLength = 500
          , removeSingletons = FALSE, keepSingletonsAbove = Inf
          , reduceToNonoverlapping = TRUE
          , useMulticore = FALSE, nrCores = NULL) {
  if (length(object) == 0) return(IRanges())
  df <- .paraclu(first(object), second(object))
  IRanges(df$start, df$end, min_d = df$min_d, max_d = df$max_d)
})

#' @rdname paraclu

setMethod("paraclu", "CTSS",
  function( object
          , minStability = 1, maxLength = 500
          , removeSingletons = FALSE, keepSingletonsAbove = Inf
          , reduceToNonoverlapping = TRUE
          , useMulticore = FALSE, nrCores = NULL) {
  # Sort and remove null ranges
  object <- sort(object)
  object <- object[score(object) != 0]
  # Split by chromosome and strand
  f <- list(seqnames(object), strand(object))
  object.list <- split(object, f, drop = TRUE)
  # Run Paraclu on each element
  object.listofpairs <- lapply(object.list, \(x) Pairs(pos(x), score(x)))
  result.list <- bplapply( object.listofpairs, paraclu
                         , BPPARAM = CAGEr_Multicore(useMulticore, nrCores))
  # Stich the results back as a GRanges object
  table <- expand.grid(seqnames = levels(f[[1]]), strand = levels(f[[2]]))
  rownames(table) <- paste(table$seqnames, table$strand, sep='.')
  clusters <- sapply(names(result.list), \(name)
    GRanges( seqnames = table[name, "seqnames"]
           , ranges   = result.list[[name]]
           , strand   = table[name, "strand"]
           , seqinfo = seqinfo(object))
  ) |> GRangesList() |> unlist()
  # Filter by stability and length
  clusters <- clusters[(clusters$max_d  >= (minStability * clusters$min_d)) &
                       (width(clusters) <= maxLength)]
  # Compute score and dominant CTSs, and remove singletons as wanted.
  clusters <-
    .ctss_summary_for_clusters( object, clusters
                              , removeSingletons = removeSingletons
                              , keepSingletonsAbove = keepSingletonsAbove)
  # Reduce to non-overlapping as wanted
  if(reduceToNonoverlapping == TRUE){
    o <- findOverlaps(clusters, drop.self = TRUE, type = "within")
    clusters <- clusters[-queryHits(o)]
  }
  # Rename mcols
  clusters$min_density <- clusters$min_d
  clusters$max_density <- clusters$max_d
  clusters$min_d <- clusters$max_d <- NULL
  # Return the clusters
  names (clusters) <- seq_along(clusters)
  as(clusters, "TagClusters")
})

setMethod("paraclu", "GRanges",
  function( object
          , minStability = 1, maxLength = 500
          , removeSingletons = FALSE, keepSingletonsAbove = Inf
          , reduceToNonoverlapping = TRUE
          , useMulticore = FALSE, nrCores = NULL) {
  stopifnot(all(width(object) == 1))
  paraclu( as(object, "CTSS"), minStability = minStability, maxLength = maxLength
         , removeSingletons = removeSingletons, keepSingletonsAbove = keepSingletonsAbove
         , reduceToNonoverlapping = reduceToNonoverlapping
         , useMulticore = useMulticore, nrCores = nrCores)
})

setMethod("paraclu", "SummarizedExperiment",
  function( object
          , minStability = 1, maxLength = 500
          , removeSingletons = FALSE, keepSingletonsAbove = Inf
          , reduceToNonoverlapping = TRUE
          , useMulticore = FALSE, nrCores = NULL) {
  
  tag.cluster.list <- GRangesList()
  for(s in colnames(object)) {
    message("\t-> ", s)
    gr <- rowRanges(object)
    score(gr) <- assays(object)[["normalizedTpmMatrix"]][[s]]
    tag.cluster.list[[s]] <-
      paraclu( gr
             , minStability = minStability, maxLength = maxLength
             , removeSingletons = removeSingletons
             , keepSingletonsAbove = keepSingletonsAbove
             , reduceToNonoverlapping = reduceToNonoverlapping
             , useMulticore = useMulticore, nrCores = nrCores)
  }
  tag.cluster.list
})
