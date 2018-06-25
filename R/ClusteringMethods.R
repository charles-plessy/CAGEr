#' @include CTSS.R Multicore.R

#' @name clusterCTSS
#' 
#' @title Cluster CTSSs into tag clusters
#' 
#' @description Clusters individual CAGE transcription start sites (CTSSs) along
#' the genome into tag clusters (TCs) using specified \emph{ab initio} method, or assigns
#' them to predefined genomic regions.
#' 
#' @param object A \code{\link{CAGEr}} object.
#' 
#' @param threshold,nrPassThreshold Ignore CTSSs with signal \code{< threshold}
#'        in \code{< nrPassThreshold} experiments.
#' 
#' @param thresholdIsTpm Logical indicating if \code{threshold} is expressed in
#'        raw tag counts (FALSE) or normalized signal (TRUE).
#' 
#' @param method Method to be used for clustering: \code{"distclu"},
#'        \code{"paraclu"} or \code{"custom"}.  See Details.
#' 
#' @param maxDist Maximal distance between two neighbouring CTSSs for them to be part of the
#'        same cluster.  Used only when \code{method = "distclu"}, otherwise ignored.
#' 
#' @param removeSingletons Logical indicating if tag clusters containing only
#'        one CTSS be removed.  Ignored when \code{method = "custom"}.
#' 
#' @param keepSingletonsAbove Controls which singleton tag clusters will be removed.
#'        When \code{removeSingletons = TRUE}, only singletons with signal
#'        \code{< keepSingletonsAbove} will be removed.  Useful to prevent removing highly
#'        supported singleton tag clusters.  Default value \code{Inf} results in removing all
#'        singleton TCs when \code{removeSingletons = TRUE}.  Ignored when
#'        \code{removeSingletons = FALSE} or \code{method = "custom"}.
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
#' @param customClusters Genomic coordinates of predefined regions to be used to
#'        segment the CTSSs.  The format is either a \code{\link{GRanges}}
#'        object or a \code{data.frame} with the following columns: \code{chr}
#'        (chromosome name), \code{start} (0-based start coordinate), \code{end}
#'        (end coordinate), \code{strand} (either \code{"+"}, or \code{"-"}).
#'        Used only when \code{method = "custom"}.
#' 
#' @param useMulticore Logical, should multicore be used.  \code{useMulticore = TRUE}
#'        has no effect on non-Unix-like platforms.
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
#' \href{http://cbrc3.cbrc.jp/~martin/paraclu/}{http://cbrc3.cbrc.jp/~martin/paraclu/}).
#' Since Paraclu finds #' clusters within clusters (unlike distclu), additional
#' parameters (\code{removeSingletons}, \code{keepSingletonsAbove},
#' \code{minStability}, \code{maxLength} and \code{reduceToNonoverlapping}) can
#' be specified to simplify the output by discarding too small (singletons) or
#' too big clusters, and to reduce the clusters to a final set of non-overlapping
#' clusters.  Clustering is done for every CAGE dataset within the CAGEr object
#' separately, resulting in a different set of tag clusters for every CAGE dataset. TCs from
#' different datasets can further be aggregated into a single referent set of consensus
#' clusters by calling the \code{\link{aggregateTagClusters}} function.
#' 
#' @return The slots \code{clusteringMethod}, \code{filteredCTSSidx} and \code{tagClusters} of
#' the provided \code{\link{CAGEset}} object will be occupied by the information on method used
#' for clustering, CTSSs included in the clusters and list of tag clusters per CAGE experiment,
#' respectively.  To retrieve tag clusters for individual CAGE dataset use
#' \code{\link{tagClusters}} function.
#' 
#' In \code{\link{CAGEexp}} objects, the results will be stored as a
#' \code{\link{GRangesList}} of \code{\link{TagClusters}} objects in the metadata
#' slot \code{tagClusters}.  The \code{TagClusters} object will contain a
#' \code{filteredCTSSidx} column if appropriate.  The clustering method name
#' is saved in the metadata slot of the \code{GRangesList}.
#' 
#' @references Frith \emph{et al.} (2007) A code for transcription initiation in mammalian
#' genomes, \emph{Genome Research} \bold{18}(1):1-12,
#' (\href{http://www.cbrc.jp/paraclu/}{http://www.cbrc.jp/paraclu/}).
#' 
#' @author Vanja Haberle
#' 
#' @seealso \code{\link{tagClusters}}, \code{\link{aggregateTagClusters}} and
#' \code{\link{CTSSclusteringMethod}}.
#' 
#' @family CAGEr object modifiers
#' @family CAGEr clusters functions
#' 
#' @examples
#' head(tagClusters(exampleCAGEset, "sample1"))
#' clusterCTSS( object = exampleCAGEset, threshold = 50, thresholdIsTpm = TRUE
#'            , nrPassThreshold = 1, method = "distclu", maxDist = 20
#'            , removeSingletons = TRUE, keepSingletonsAbove = 100)
#' head(tagClusters(exampleCAGEset, "sample1"))
#' 
#' clusterCTSS( exampleCAGEexp, threshold = 50, thresholdIsTpm = TRUE
#'            , nrPassThreshold = 1, method = "distclu", maxDist = 20
#'            , removeSingletons = TRUE, keepSingletonsAbove = 100)
#' tagClustersGR(exampleCAGEexp, "Zf.30p.dome")
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

setMethod( "clusterCTSS", "CAGEset"
         , function( object, threshold, nrPassThreshold, thresholdIsTpm, method
                   , maxDist, removeSingletons, keepSingletonsAbove, minStability
                   , maxLength, reduceToNonoverlapping, customClusters
                   , useMulticore, nrCores) {
  
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
	
	getCTSSdataSE <- function() {
	if(threshold > 0){
		nr.pass.threshold <- apply(data, 1, function(x) {sum(x >= threshold)})
		idx <<- nr.pass.threshold >= min(nrPassThreshold, length(sample.labels))
		data <<- CTSStagCountSE(object)[idx,]
	}else{
		data <<- CTSStagCountSE(object)
		idx <<- rep(TRUE, nrow(data))
	}}
	
	getCTSSdatadf <- function() {
	if(threshold > 0){
		nr.pass.threshold <- apply(data, 1, function(x) {sum(x >= threshold)})
		idx <<- nr.pass.threshold >= min(nrPassThreshold, length(sample.labels))
		data <<- cbind(CTSScoordinates(object)[idx,], object@normalizedTpmMatrix[idx,,drop=F])
	}else{
		data <<- cbind(CTSScoordinates(object), object@normalizedTpmMatrix)
		idx <<- rep(TRUE, nrow(data))
	}}
		
	message("Clustering...")
	method <- match.arg(method)
	if(method == "distclu"){
	  getCTSSdataSE()
		ctss.cluster.list <- .distclu(se = data, max.dist = maxDist, removeSingletons = removeSingletons, keepSingletonsAbove = keepSingletonsAbove, useMulticore = useMulticore, nrCores = nrCores)
	}else if (method == "paraclu"){
	  getCTSSdatadf()
		ctss.cluster.list <- .paraclu(data = data, sample.labels = sample.labels, minStability = minStability, maxLength = maxLength, removeSingletons = removeSingletons, keepSingletonsAbove = keepSingletonsAbove, reduceToNonoverlapping = reduceToNonoverlapping, useMulticore = useMulticore, nrCores = nrCores)
	}else if(method == "custom"){
	  getCTSSdatadf()
		if(length(customClusters)==0){
			stop("'customClusters' must be given when method = \"custom\"")
		}
		ctss.cluster.list <- .predefined.clusters(data = data, sample.labels = sample.labels, custom.clusters = customClusters, useMulticore = useMulticore, nrCores = nrCores)
	}
	
	object@filteredCTSSidx <- idx
	object@clusteringMethod <- method
	if (method == "distclu") {
	  ctss.cluster.list <- lapply(ctss.cluster.list, TCgranges2dataframe)
	}
	object@tagClusters <- ctss.cluster.list
	assign(objName, object, envir = parent.frame())
	invisible(1)
})

#' @rdname clusterCTSS

setMethod( "clusterCTSS", "CAGEexp"
         , function( object, threshold, nrPassThreshold, thresholdIsTpm, method, maxDist
                   , removeSingletons, keepSingletonsAbove, minStability, maxLength
                   , reduceToNonoverlapping, customClusters, useMulticore, nrCores) {

  objName <- deparse(substitute(object))
  assay <- ifelse(thresholdIsTpm, "normalizedTpmMatrix", "counts")
  data <- CTSStagCountSE(object)

  if (! "normalizedTpmMatrix" %in% names(assays(data)))
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
    ctss.cluster.list <- .paraclu( data = data[decode(filteredCTSSidx(object)),]
                                 , sample.labels = sampleLabels(object)
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
  metadata(object)$tagClusters <- ctss.cluster.list
  assign(objName, object, envir = parent.frame())
  invisible(1)
})

#' @name .clusterAggregateAndSum
#' @rdname clusterAggregateAndSum
#' 
#' @param clusters Clusters to be aggregated.  \code{data.frame}, or
#' \code{ConsensusClusters}, which will be coerced to \code{data.frame}.
#' 
#' @param key Name of the column containing the factor used to aggregate
#' the clusters.
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
