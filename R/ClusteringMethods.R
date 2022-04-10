#' @include CTSS.R Multicore.R

#' @name clusterCTSS
#' 
#' @title Cluster CTSSs into tag clusters
#' 
#' @description Clusters individual CAGE transcription start sites (CTSSs) along
#' the genome into tag clusters (TCs) using specified _ab initio_ method, or
#' assigns them to predefined genomic regions.
#' 
#' @param object A [`CAGEr`] object.
#' 
#' @param threshold,nrPassThreshold Ignore CTSSs with signal `< threshold`
#'        in `< nrPassThreshold` experiments.
#' 
#' @param thresholdIsTpm Logical indicating if `threshold` is expressed in
#'        raw tag counts (`FALSE`) or normalized signal (`TRUE`).
#' 
#' @param method Method to be used for clustering: `"distclu"`,
#'        `"paraclu"` or `"custom"`.  See Details.
#' 
#' @param maxDist Maximal distance between two neighbouring CTSSs for them to be
#'        part of the same cluster.  Used only when `method = "distclu"`,
#'        otherwise ignored.
#' 
#' @param removeSingletons Logical indicating if tag clusters containing only
#'        one CTSS be removed.  Ignored when `method = "custom"`.
#' 
#' @param keepSingletonsAbove Controls which singleton tag clusters will be
#'        removed.  When `removeSingletons = TRUE`, only singletons with signal
#'        `< keepSingletonsAbove` will be removed.  Useful to prevent removing
#'        highly supported singleton tag clusters.  Default value `Inf` results
#'        in removing all singleton TCs when `removeSingletons = TRUE`.  Ignored
#'        when `method = "custom"`.
#' 
#' @param minStability Minimal stability of the cluster, where stability is
#'        defined as ratio between maximal and minimal density value for which
#'        this cluster is maximal scoring.  For definition of stability refer to
#'        Frith _et al._, Genome Research, 2007.  Clusters with stability
#'        `< minStability` will be discarded.  Used only when `method = "paraclu"`.
#' 
#' @param maxLength Maximal length of cluster in base-pairs.  Clusters with length
#'        `> maxLength` will be discarded.  Ignored when `method = "custom"`.
#' 
#' @param reduceToNonoverlapping Logical, should smaller clusters contained
#'        within bigger cluster be removed to make a final set of tag clusters
#'        non-overlapping. Used only `method = "paraclu"`.
#' 
#' @param customClusters Genomic coordinates of predefined regions to be used to
#'        segment the CTSSs.  The format is either a [`GRanges`] object or a
#'        [`data.frame`] with the following columns: `chr` (chromosome name),
#'        `start` (0-based start coordinate), `end` (end coordinate), `strand`
#'        (either `"+"`, or `"-"`).  Used only when `method = "custom"`.
#' 
#' @param useMulticore Logical, should multicore be used.  `useMulticore = TRUE`
#'        has no effect on non-Unix-like platforms.
#' 
#' @param nrCores Number of cores to use when `useMulticore = TRUE`.  Default
#'        value `NULL` uses all detected cores.
#' 
#' @details The `"distclu"` method is an implementation of simple distance-based
#' clustering of data attached to sequences, where two neighbouring TSSs are
#' joined together if they are closer than some specified distance (see
#' [`distclu-functions`] for implementation details.
#' 
#' `"paraclu"` is an implementation of Paraclu algorithm for parametric
#' clustering of data attached to sequences (Frith _et al._, Genome Research,
#' 2007).  Since Paraclu finds clusters within clusters (unlike distclu),
#' additional parameters (`removeSingletons`, `keepSingletonsAbove`,
#' `minStability`, `maxLength` and `reduceToNonoverlapping`) can be specified to
#' simplify the output by discarding too small (singletons) or too big clusters,
#' and to reduce the clusters to a final set of non-overlapping clusters.
#' 
#' Clustering is done for every CAGE dataset within the CAGEr object separately,
#' resulting in a different set of tag clusters for every CAGE dataset. TCs from
#' different datasets can further be aggregated into a single referent set of
#' consensus clusters by calling the [`aggregateTagClusters`] function.
#' 
#' @return Returns the [`CAGEexp`] object, in which, the results will be stored as a `GRangesList` of
#' [`TagClusters`] objects in the metadata slot `tagClusters`.  The
#' `TagClusters` objects will contain a `filteredCTSSidx` column if appropriate.
#' The clustering method name is saved in the metadata slot of the `GRangesList`.
#' 
#' @references Frith _et al._ (2007) A code for transcription initiation in
#' mammalian genomes, _Genome Research_ **18**(1):1-12,
#' (\href{http://www.cbrc.jp/paraclu/}{http://www.cbrc.jp/paraclu/}).
#' 
#' @author Vanja Haberle
#' 
#' @seealso [`tagClustersGR`], [`aggregateTagClusters`] and [`CTSSclusteringMethod`].
#' 
#' @family CAGEr object modifiers
#' @family CAGEr clusters functions
#' 
#' @examples
#' 
#' # Using 'distclu', notice argument 'maxDist'
#' ce <- clusterCTSS( exampleCAGEexp, threshold = 50, thresholdIsTpm = TRUE
#'            , nrPassThreshold = 1, method = "distclu", maxDist = 20
#'            , removeSingletons = TRUE, keepSingletonsAbove = 100)
#' tagClustersGR(ce, "Zf.30p.dome")
#' 
#' # Using 'paraclu', notice arguments 'maxLength' and 'minStability'
#' ce <- clusterCTSS( exampleCAGEexp, threshold = 50, thresholdIsTpm = TRUE
#'            , nrPassThreshold = 1, method = "paraclu"
#'            , removeSingletons = TRUE, keepSingletonsAbove = 100
#'            , maxLength = 500, minStability = 1
#'            , reduceToNonoverlapping = TRUE)
#' tagClustersGR(ce, "Zf.30p.dome")
#' 
#' @export

setGeneric( "clusterCTSS"
          , function( object
                    , threshold = 1, nrPassThreshold = 1, thresholdIsTpm = TRUE
                    , method = c("distclu", "paraclu", "custom"), maxDist = 20
                    , removeSingletons = FALSE, keepSingletonsAbove = Inf
                    , minStability = 1, maxLength = 500
                    , reduceToNonoverlapping = TRUE, customClusters = NULL
                    , useMulticore = FALSE, nrCores = NULL)
              standardGeneric("clusterCTSS"))

#' @rdname clusterCTSS

setMethod( "clusterCTSS", "CAGEexp"
         , function( object, threshold, nrPassThreshold, thresholdIsTpm, method, maxDist
                   , removeSingletons, keepSingletonsAbove, minStability, maxLength
                   , reduceToNonoverlapping, customClusters, useMulticore, nrCores) {

  assay <- ifelse(isTRUE(thresholdIsTpm), "normalizedTpmMatrix", "counts")
  data <- CTSStagCountSE(object)

  if (! "normalizedTpmMatrix" %in% assayNames(data))
    stop( "Could not find normalized CAGE signal values, see ?normalizeTagCount.\n"
        , "clusterCTSS() needs normalized values to create its output tables, that "
        , "include TPM expression columns.")

  message("\nFiltering out CTSSs below threshold...")
  filteredCTSSidx(object) <-
    .filterCtss(data, threshold = threshold
               , nrPassThreshold = nrPassThreshold, thresholdIsTpm = thresholdIsTpm)
	
	message("Clustering...")
	method <- match.arg(method)

  if (method == "distclu") {
    ctss.cluster.list <- .distclu( se = data[decode(filteredCTSSidx(object)),]
                                 , max.dist = maxDist, removeSingletons = removeSingletons
                                 , keepSingletonsAbove = keepSingletonsAbove
                                 , useMulticore = useMulticore, nrCores = nrCores)
  } else if (method == "paraclu") {
    ctss.cluster.list <- .paraclu( se = data[decode(filteredCTSSidx(object)),]
                                 , minStability = minStability, maxLength = maxLength
                                 , removeSingletons = removeSingletons
                                 , keepSingletonsAbove = keepSingletonsAbove
                                 , reduceToNonoverlapping = reduceToNonoverlapping
                                 , useMulticore = useMulticore, nrCores = nrCores)
  } else if(method == "custom") {
    if(is.null(customClusters))
    	stop(sQuote("customClusters"), " must be given when method = ", sQuote("custom"), ".")
    ctss.cluster.list <- .predefined.clusters( data = data[decode(filteredCTSSidx(object)),]
                                             , sample.labels = sampleLabels(object)
                                             , custom.clusters = customClusters
                                             , useMulticore = useMulticore, nrCores = nrCores)
  }
  
  CTSSclusteringMethod(ctss.cluster.list) <- method
  seqlevels(ctss.cluster.list) <- seqlevels(CTSStagCountSE(object))
  seqinfo(ctss.cluster.list)   <- seqinfo(CTSStagCountSE(object))
  metadata(object)$tagClusters <- ctss.cluster.list
  object
})

#' @name .clusterAggregateAndSum
#' @rdname clusterAggregateAndSum
#' 
#' @param clusters Clusters to be aggregated.  `data.frame`, or
#' `GRanges`, which will be coerced to `data.frame`.
#' 
#' @param key Name of the column containing the factor used to aggregate
#' the clusters.
#' 
#' @title Aggregate identical clusters and sum their scores.
#' 
#' @description Private function using `data.table` objects to preform grouping
#' operations at a high performance.  These functions use _non-standard
#' evaluation_ in a context that raises warnings in `R CMD check`.  By
#' separating these functions from the rest of the code, I hope to make the
#' workarounds easier to manage.

setGeneric(".clusterAggregateAndSum", function (clusters, key) standardGeneric(".clusterAggregateAndSum"))

#' @rdname clusterAggregateAndSum
#' @importFrom data.table setkeyv setnames

setMethod(".clusterAggregateAndSum", "data.table", function (clusters, key) {
  setkeyv(clusters, key)
  chr <- min <- max <- strand <- tpm <- NULL
	clusters <- clusters[ , list( chr[1]
	                            , min(start)
	                            , max(end)
	                            , strand[1]
	                            , sum(tpm))
	                      , by = key]
	setnames(clusters, c(key, "chr", "start", "end", "strand", "tpm"))
	setkeyv(clusters, key)
})

#' @rdname clusterAggregateAndSum
#' @importFrom data.table data.table

setMethod(".clusterAggregateAndSum", "data.frame", function (clusters, key) {
  as.data.frame(.clusterAggregateAndSum(data.table(clusters), key))
})

#' @rdname clusterAggregateAndSum

setMethod(".clusterAggregateAndSum", "GRanges", function (clusters, key) {
  CCdataframe2granges(.clusterAggregateAndSum(CCgranges2dataframe(clusters), key))
})


#' @rdname byCtss
#' 
#' @title Apply functions to identical CTSSes.
#' 
#' @param ctssDT A \code{\link{data.table}} representing CTSSes.
#' @param colName The name of the column on which to apply the function.
#' @param fun The function to apply.
#' 
#' @description \code{.byCTSS} is a private function using  \code{data.table} objects
#' to preform grouping operations at a high performance.  These functions use
#' \emph{non-standard evaluation} in a context that raises warnings in \code{R CMD check}.
#' By separating these functions from the rest of the code, I hope to make the workarounds
#' easier to manage.
#' 
#' @examples
#' ctssDT <- data.table::data.table(
#'   chr       = c("chr1", "chr1", "chr1", "chr2"),
#'   pos       = c(1     , 1     , 2     , 1     ),
#'   strand    = c("+"   , "+"   , "-"   , "-"   ),
#'   tag_count = c(1     , 1     , 1     , 1     ))
#' ctssDT
#' CAGEr:::.byCtss(ctssDT, "tag_count", sum)

setGeneric( ".byCtss"
          , function (ctssDT, colName, fun) standardGeneric(".byCtss"))

#' @rdname byCtss

setMethod(".byCtss", "data.table", function (ctssDT, colName, fun) {
  if (! all(c("chr", "pos", "strand") %in% colnames(ctssDT))) stop("These are not CTSSes.")
  chr <- pos <- strand <- .SD <- NULL
  ctssDT[ , fun(.SD[[1]])
        , by = list(chr, pos, strand)
        , .SDcols = colName]
})
