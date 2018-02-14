#' @name CTSS-class
#' @aliases CTSS
#' 
#' @title Classes for type safety.
#' 
#' @description The `CTSS` class represents CAGE transcription start sites
#' (CTSS) at single-nucleotide resolution, using `GPos` as base class.  It is
#' used internally by _CAGEr_ for type safety.
#' 
#' @details The `genomeName` element of the `metadata` slot is used to store the
#' name of the _BSgenome_ package used when constructing the `CAGEr` object.
#' Be careful that calling the `CTSS` constructor twice in a row will erase
#' the fist metadata (because the `GPos` constructor does so).
#' 
#' @rdname CTSS-class
#' 
#' @import methods
#' @importFrom GenomicRanges GRanges
#' 
#' @author Charles Plessy
#' 
#' @examples 
#' ctss <- CAGEr:::.CTSS( CTSScoordinatesGR(exampleCAGEexp)
#'                      , bsgenomeName = genomeName(exampleCAGEexp))
#' genomeName(ctss)
#' 
#' @export

# (See <https://github.com/Bioconductor/Contributions/issues/261#issuecomment-277479436>.)

.CTSS <- setClass("CTSS", contains = "GPos")

setMethod("initialize", "CTSS", function(.Object, ..., bsgenomeName = NULL) {
  .Object <- callNextMethod(.Object, ...)
  if (! is.null(bsgenomeName))
    metadata(.Object)$genomeName <- bsgenomeName
  .Object
})

#' @name CTSS.chr-class
#' 
#' @aliases CTSS.chr
#' 
#' @details The `CTSS.chr` class represents a [CTSS] object that is guaranteed
#' to be only on a single chromosome.  It is used internally by _CAGEr_ for type
#' safety.
#' 
#' @rdname CTSS-class
#' 
#' @import methods
#' @importFrom GenomicRanges GPos
#' 
#' @examples 
#' ctss.chr <- as(CTSScoordinatesGR(exampleCAGEexp), "CTSS.chr")

.CTSS.chr <-
  setClass( "CTSS.chr"
          , contains = "CTSS"
          , validity =
  function(object) 
    if (length(unique(seqnames(object))) != 1)
      return("Mutiple sequnames found: CTSS.chr objects should be only on a single chromosome.")
)


#' @name ConsensusClusters
#' 
#' @rdname CTSS-class
#' 
#' @details The \code{ConsensusClusters} class represents consensus clusters.
#' It is used internally by CAGEr for type safety.

.ConsensusClusters <- setClass( "ConsensusClusters", contains = "GRanges")

#' @name TagClusters
#' 
#' @rdname CTSS-class
#' 
#' @details The \code{TagClusters} class represents tag clusters.
#' It is used internally by CAGEr for type safety.

.TagClusters <- setClass( "TagClusters", contains = "GRanges")
