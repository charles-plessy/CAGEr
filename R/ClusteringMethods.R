#' @include CTSS.R Multicore.R

#' @name clusterCTSS
#' 
#' @title Clustering CTSSs into tag clusters
#' 
#' @description Clusters individual CAGE transcription start sites (CTSSs) along
#' the genome into tag clusters (TCs) using specified \emph{ab initio} method, or assigns
#' them to predefined genomic regions.
#' 
#' @param object A \code{\link{CAGEr}} object.
#' 
#' @param threshold,nrPassThreshold Only CTSSs with signal \code{>= threshold} in
#'        \code{>= nrPassThreshold} experiments will be used for clustering and will
#'        contribute towards total signal of the cluster.
#' 
#' @param thresholdIsTpm Logical, is threshold raw tag count value (FALSE) or
#'        normalized signal (TRUE).
#' 
#' @param method Method to be used for clustering. Can be one of the \code{"distclu"},
#'        \code{"paraclu"} or \code{"custom"}.  See Details.
#' 
#' @param maxDist Maximal distance between two neighbouring CTSSs for them to be part of the
#'        same cluster.  Used only when \code{method = "distclu"}, otherwise ignored.
#' 
#' @param removeSingletons Logical, should tag clusters containing only one CTSS be removed.
#'        Ignored when \code{method = "custom"}.
#' 
#' @param keepSingletonsAbove Controls which singleton tag clusters will be removed.
#'        When \code{removeSingletons = TRUE}, only singletons with signal
#'        \code{< keepSingletonsAbove} will be removed.  Useful to prevent removing highly
#'        supported singleton tag clusters.  Default value \code{Inf} results in removing all
#'         singleton TCs when \code{removeSingletons = TRUE}.  Ignored when
#'         \code{removeSingletons = FALSE} or \code{method = "custom"}.
#' 
#' @param minStability Minimal stability of the cluster, where stability is defined as ratio
#'        between maximal and minimal density value for which this cluster is maximal scoring.
#'        For definition of stability refer to Frith \emph{et al}., Genome Research, 2007.
#'        Clusters with stability \code{< minStability} will be discarded.  Used only when
#'        \code{method = "paraclu"}, otherwise ignored.
#' 
#' @param maxLength Maximal length of cluster in base-pairs.  Clusters with length
#'        \code{> maxLength} will be discarded.  Ignored when \code{method = "custom"}.
#' 
#' @param reduceToNonoverlapping Logical, should smaller clusters contained within bigger
#'        cluster be removed to make a final set of tag clusters non-overlapping. Used only
#'        when \code{method = "paraclu"}.  See Details.
#' 
#' @param customClusters Genomic coordinates of predefined regions to be used to segment the
#'        CTSSs.  It has to be a \code{data.frame} with following columns: \code{chr}
#'        (chromosome name), \code{start} (0-based start coordinate), \code{end}
#'        (end coordinate), \code{strand} (either \code{"+"}, or \code{"-"}).  Used only when
#'        \code{method = "custom"}.
#' 
#' @param useMulticore Logical, should multicore be used.  \code{useMulticore = TRUE} is
#'        supported only on Unix-like platforms.
#' 
#' @param nrCores Number of cores to use when \code{useMulticore = TRUE}.  Default value
#'        \code{NULL} uses all detected cores.
#' 
#' @details Two \emph{ab initio} methods for clustering TSSs along the genome are supported:
#' \code{"distclu"} and \code{"paraclu"}.  \code{"distclu"} is an implementation of simple
#' distance-based clustering of data attached to sequences, where two neighbouring TSSs are
#' joined together if they are closer than some specified distance (see
#' \code{\link{distclu-functions}} for implementation details.  \code{"paraclu"} is an
#' implementation of Paraclu algorithm for parametric clustering of data attached to
#' sequences developed by M. Frith (Frith \emph{et al.}, Genome Research, 2007,
#' \href{http://www.cbrc.jp/paraclu/}{http://www.cbrc.jp/paraclu/}).  Since Paraclu finds
#' clusters within clusters (unlike distclu), additional parameters (\code{removeSingletons},
#' \code{keepSingletonsAbove}, \code{minStability}, \code{maxLength} and
#' \code{reduceToNonoverlapping}) can be specified to simplify the output by discarding too
#' small (singletons) or too big clusters, and to reduce the clusters to a final set of
#' non-overlapping clusters.  Clustering is done for every CAGE dataset within CAGEset object
#' separatelly, resulting in a different set of tag clusters for every CAGE dataset. TCs from
#' different datasets can further be aggregated into a single referent set of consensus
#' clusters by calling \code{\link{aggregateTagClusters}} function.
#' 
#' @return The slots \code{clusteringMethod}, \code{filteredCTSSidx} and \code{tagClusters} of
#' the provided \code{\link{CAGEset}} object will be occupied by the information on method used
#' for clustering, CTSSs included in the clusters and list of tag clusters per CAGE experiment,
#' respectively.  To retrieve tag clusters for individual CAGE dataset use
#' \code{\link{tagClusters}} function.
#' 
#' @references Frith \emph{et al.} (2007) A code for transcription initiation in mammalian
#' genomes, \emph{Genome Research} \bold{18}(1):1-12,
#' (\href{http://www.cbrc.jp/paraclu/}{http://www.cbrc.jp/paraclu/}).
#' 
#' @author Vanja Haberle
#' 
#' @seealso \code{\link{tagClusters}}, \code{\link{aggregateTagClusters}}
#' 
#' @family CAGEr object modifiers
#' @family CAGEr clusters functions
#' 
#' @examples
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' head(tagClusters(exampleCAGEset, "sample1"))
#' clusterCTSS( object = exampleCAGEset, threshold = 50, thresholdIsTpm = TRUE
#'            , nrPassThreshold = 1, method = "distclu", maxDist = 20
#'            , removeSingletons = TRUE, keepSingletonsAbove = 100)
#' head(tagClusters(exampleCAGEset, "sample1"))
#' 
#' ce <- readRDS(system.file(package = "CAGEr", "extdata/CAGEexp.rds"))
#' normalizeTagCount(ce)
#' clusterCTSS( object = ce, threshold = 50, thresholdIsTpm = TRUE
#'            , nrPassThreshold = 1, method = "distclu", maxDist = 20
#'            , removeSingletons = TRUE, keepSingletonsAbove = 100)
#' head(tagClusters(ce, "Zf.30p.dome"))
#' 
#' @export

setGeneric( "clusterCTSS"
          , function( object
                    , threshold = 1, nrPassThreshold = 1, thresholdIsTpm = TRUE
                    , method = "distclu", maxDist = 20
                    , removeSingletons = FALSE, keepSingletonsAbove = Inf
                    , minStability = 1, maxLength = 500
                    , reduceToNonoverlapping = TRUE, customClusters = NULL
                    , useMulticore = FALSE, nrCores = NULL)
              standardGeneric("clusterCTSS"))

#' @rdname clusterCTSS

setMethod("clusterCTSS",
signature(object = "CAGEset"),
function (object, threshold, nrPassThreshold, thresholdIsTpm, method, maxDist, removeSingletons, keepSingletonsAbove, minStability, maxLength, reduceToNonoverlapping, customClusters, useMulticore, nrCores){
  
  useMulticore <- .checkMulticore(useMulticore)
		
	objName <- deparse(substitute(object))
	sample.labels <- sampleLabels(object)
	
	message("\nFiltering CTSSs below threshold...")
	if(thresholdIsTpm){
	  if (identical(object@normalizedTpmMatrix, data.frame()))
	    stop("Could not find normalized CAGE signal values, see ?normalizeTagCount.")
		data <- object@normalizedTpmMatrix
	}else{
		if (identical(object@tagCountMatrix, data.frame()))
	    stop("Could not find CTSS tag counts, see ?getCTSS.")
	  data <- object@tagCountMatrix
	}
	
	if (identical(object@normalizedTpmMatrix, data.frame()))
	    stop("Could not find normalized CAGE signal values, see ?normalizeTagCount.\n",
	         "clusterCTSS() needs normalized values to create its output tables, that ",
	         "include TPM expression columns.")
	
	if(threshold > 0){
		nr.pass.threshold <- apply(data, 1, function(x) {sum(x >= threshold)})
		idx <- nr.pass.threshold >= min(nrPassThreshold, length(sample.labels))
		gc()
		data <- CTSStagCountSE(object)[idx,]
		gc()
	}else{
		data <- CTSStagCountSE(object)
		idx <- rep(TRUE, nrow(data))

	}
	
	message("Clustering...")
	if(method == "distclu"){
		ctss.cluster.list <- .distclu(se = data, max.dist = maxDist, removeSingletons = removeSingletons, keepSingletonsAbove = keepSingletonsAbove, useMulticore = useMulticore, nrCores = nrCores)
	}else if (method == "paraclu"){
		ctss.cluster.list <- .paraclu(data = data, sample.labels = sample.labels, minStability = minStability, maxLength = maxLength, removeSingletons = removeSingletons, keepSingletonsAbove = keepSingletonsAbove, reduceToNonoverlapping = reduceToNonoverlapping, useMulticore = useMulticore, nrCores = nrCores)
	}else if(method == "custom"){
		if(length(customClusters)==0){
			stop("'customClusters' must be given when method = \"custom\"")
		}
		ctss.cluster.list <- .predefined.clusters(data = data, sample.labels = sample.labels, custom.clusters = customClusters, useMulticore = useMulticore, nrCores = nrCores)
	}else{
		stop("'method' parameter must be one of the (\"distclu\", \"paraclu\", \"custom\")")
	}
	
	object@filteredCTSSidx <- idx
	object@clusteringMethod <- method
	object@tagClusters <- lapply(ctss.cluster.list, TCgranges2dataframe)
	assign(objName, object, envir = parent.frame())
	invisible(1)
})

#' @rdname clusterCTSS

setMethod("clusterCTSS",
signature(object = "CAGEexp"),
function (object, threshold, nrPassThreshold, thresholdIsTpm, method, maxDist, removeSingletons, keepSingletonsAbove, minStability, maxLength, reduceToNonoverlapping, customClusters, useMulticore, nrCores){
  
  useMulticore <- .checkMulticore(useMulticore)
  objName <- deparse(substitute(object))
  assay <- ifelse(thresholdIsTpm, "normalizedTpmMatrix", "counts")

  message("\nFiltering CTSSs below threshold...")

  if (! "tagCountMatrix" %in% names(experiments(object)))
    stop("Could not find CTSS tag counts, see ?getCTSS.")

  if (! "normalizedTpmMatrix" %in% names(assays(CTSStagCountSE(object))))
    stop( "Could not find normalized CAGE signal values, see ?normalizeTagCount.\n"
        , "clusterCTSS() needs normalized values to create its output tables, that "
        , "include TPM expression columns.")

  data <- CTSStagCountSE(object)

	if (threshold > 0) {
		nr.pass.threshold <- rowSums(DelayedArray(assays(data)[[assay]]) >= threshold)
		rowRanges(data)$filteredCTSSidx <-
		  Rle(nr.pass.threshold >= min(nrPassThreshold, length(sampleLabels(object))))
	} else {
	  rowRanges(data)$filteredCTSSidx <- Rle(TRUE)
	}
	
	message("Clustering...")
	
  if (method == "distclu") {
    ctss.cluster.list <- .distclu( se = data[rowRanges(data)$filteredCTSSidx,]
                                 , max.dist = maxDist, removeSingletons = removeSingletons
                                 , keepSingletonsAbove = keepSingletonsAbove
                                 , useMulticore = useMulticore, nrCores = nrCores)
  } else if (method == "paraclu") {
    ctss.cluster.list <- .paraclu( data = data[rowRanges(data)$filteredCTSSidx,]
                                 , sample.labels = sampleLabels(object)
                                 , minStability = minStability, maxLength = maxLength
                                 , removeSingletons = removeSingletons
                                 , keepSingletonsAbove = keepSingletonsAbove
                                 , reduceToNonoverlapping = reduceToNonoverlapping
                                 , useMulticore = useMulticore, nrCores = nrCores)
  } else if(method == "custom") {
    if(is.null(customClusters))
    	stop("'customClusters' must be given when method = \"custom\"")
    ctss.cluster.list <- .predefined.clusters( data = data[rowRanges(data)$filteredCTSSidx,]
                                             , sample.labels = sampleLabels(object)
                                             , custom.clusters = customClusters
                                             , useMulticore = useMulticore, nrCores = nrCores)
  } else {
    stop( sQuote("method"), " parameter must be "
        , dQuote("distclu"), ", ", dQuote("paraclu"), ", or ", dQuote("custom"), ".")
  }
  
	CTSStagCountSE(object) <- data
  metadata(object)$clusteringMethod <- method
  metadata(object)$tagClusters <- ctss.cluster.list
  assign(objName, object, envir = parent.frame())
  invisible(1)
})

#' @name .clusterAggregateAndSum
#' @rdname clusterAggregateAndSum
#' 
#' @title Aggregate identical clusters and sum their scores.
#' 
#' @description Private function using  \code{data.table} objects to preform grouping
#' operations at a high performance.  These functions use \emph{non-standard evaluation}
#' in a context that raises warnings in \code{R CMD check}.  By separating these functions
#' from the rest of the code, I hope to make the workarounds easier to manage.

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

setMethod(".clusterAggregateAndSum", "ConsensusClusters", function (clusters, key) {
  CCdataframe2granges(.clusterAggregateAndSum(CCgranges2dataframe(clusters), key))
})


#' @name .ctssAggregateAndSum
#' @rdname ctssAggregateAndSum
#' 
#' @title Aggregate identical CTSS and sum their scores.
#' 
#' @description Private function using  \code{data.table} objects to preform grouping
#' operations at a high performance.  These functions use \emph{non-standard evaluation}
#' in a context that raises warnings in \code{R CMD check}.  By separating these functions
#' from the rest of the code, I hope to make the workarounds easier to manage.

setGeneric(".ctssAggregateAndSum", function (ctss) standardGeneric(".ctssAggregateAndSum"))

#' @rdname ctssAggregateAndSum
#' @importFrom data.table setkeyv setnames
#' @example 
#' ctssDT <- data.table( chr       = c("chr1", "chr1", "chr1", "chr2")
#'                     , pos       = c(1     , 1     , 2     , 1     )
#'                     , strand    = c("+"   , "+"   , "-"   , "-"   )
#'                     , tag_count = c(1     , 1     , 1     , 1     ))
#' ctssDT
#' .ctssAggregateAndSum(ctssDT)

setMethod(".ctssAggregateAndSum", "data.table", function (ctss) {
  if (! all(c("chr", "pos", "strand") %in% colnames(ctss))) stop("These are not CTSSes.")
  if (! "tag_count" %in% colnames(ctss)) stop("Missing ", sQuote("tag_count"), " column.")
  chr <- pos <- strand <- tag_count <- NULL
  ctssDT[ , as.integer(sum(tag_count))
          , by = list(chr, pos, strand)]
})