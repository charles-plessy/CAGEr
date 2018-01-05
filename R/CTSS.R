#' @name CTSS-class
#' @aliases CTSS
#' 
#' @title Classes for type safety.
#' 
#' @description  The \code{CTSS} class represents CAGE transcription start sites (CTSS) at
#' single-nucleotide resolution, using GRanges as base class.  It is used internally by CAGEr
#' for type safety.
#' 
#' @rdname CTSS-class
#' 
#' @import methods
#' @importFrom GenomicRanges GRanges
#' 
#' @author Charles Plessy
#' 
#' @examples 
#' ctss <- CAGEr:::.CTSS(CTSScoordinatesGR(exampleCAGEexp))
#' 
#' @export

# (See <https://github.com/Bioconductor/Contributions/issues/261#issuecomment-277479436>.)

.CTSS <-
  setClass( "CTSS"
          , contains = "GRanges"
          , validity = 
  function(object)
    if (! identical(start(object), end(object)))
    return("Not a CTSS: start and end positions differ.")
)

#' @name CTSS.chr-class
#' 
#' @details The \code{CTSS.chr} class represents a \code{CTSS} object that is guaranteed to be
#' only on a single chromosome.  It is used internally by CAGEr for type safety.
#' 
#' @rdname CTSS-class
#' 
#' @import methods
#' 
#' @examples 
#' ctss.chr <- CAGEr:::.CTSS.chr(CTSScoordinatesGR(exampleCAGEexp))

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
