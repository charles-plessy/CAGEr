#' @include GetMethods.R

#' @name cumulativeCTSSdistribution
#' 
#' @title Cumulative sums of CAGE counts along genomic regions
#' 
#' @description Calculates the cumulative sum of normalised CAGE counts along each tag
#' cluster or consensus cluster in every sample within a CAGEr object.
#' 
#' @param object A [`CAGEr`] object
#' 
#' @param clusters `tagClusters` or `consensusClusters`.
#'        
#' @param useMulticore Logical, should multicore be used.
#'        `useMulticore = TRUE` has no effect on non-Unix-like platforms.
#' 
#' @param nrCores Number of cores to use when `useMulticore = TRUE`
#'        (set to `NULL` to use all detected cores).
#' 
#' @return In `CAGEexp` objects, cumulative sums for the _tag clusters_ are
#' stored in the metadata slot using the `RleList` class.  For _consensus
#' clusters_, they are stored in _assays_ of the `consensusClusters` experiment
#' slot of the `CAGEexp` object.
#' 
#' @author Vanja Haberle
#' @author Charles Plessy
#' 
#' @family CAGEr object modifiers
#' @family CAGEr clusters functions
#' 
#' @examples
#' cumulativeCTSSdistribution(exampleCAGEexp, clusters = "tagClusters")
#' CTSScumulativesTagClusters(exampleCAGEexp)[[1]][1:6]
#' cumulativeCTSSdistribution(exampleCAGEexp, clusters = "consensusClusters")
#' CTSScumulativesCC(exampleCAGEexp)[[1]][1:6]
#'            
#' @importFrom IRanges RleList
#' @export

setGeneric( "cumulativeCTSSdistribution"
          , function ( object, clusters = c("tagClusters", "consensusClusters")
                     , useMulticore = FALSE, nrCores = NULL)
              standardGeneric("cumulativeCTSSdistribution"))

#' @rdname cumulativeCTSSdistribution

setMethod( "cumulativeCTSSdistribution", "CAGEexp"
         , function (object, clusters, useMulticore, nrCores) {
  message("\nCalculating cumulative sum of CAGE signal along clusters...")
  clusters <- match.arg(clusters)
  getClusters <- switch( clusters
                       , tagClusters       = tagClustersGR
                       , consensusClusters = consensusClustersGR)
  setClusters <- switch(clusters
                       , tagClusters       = `CTSScumulativesTagClusters<-`
                       , consensusClusters = `CTSScumulativesCC<-`)
  samples.cumsum.list <- lapply(sampleList(object), function(s) {
		message("\t-> ", s)
    ctss <- CTSSnormalizedTpmGR(object, s)
    if (!is.null(ctss$filteredCTSSidx))
      ctss <- ctss[ctss$filteredCTSSidx]
    .getCumsum( ctss         = ctss
              , clusters     = getClusters(object, s)
              , useMulticore = useMulticore, nrCores = nrCores)
  })
  samples.cumsum.list <- lapply(samples.cumsum.list, RleList)
  setClusters(object, samples.cumsum.list)
})
