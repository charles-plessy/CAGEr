#' @name CTSS-class
#' 
#' @title CAGE Transcription Start Site class
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
#' ce <- readRDS(system.file(package = "CAGEr", "extdata/CAGEexp.rds"))
#' ctss <- CAGEr:::.CTSS(CTSScoordinatesGR(ce))
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
#' ce <- readRDS(system.file(package = "CAGEr", "extdata/CAGEexp.rds"))
#' ctss.chr <- CAGEr:::.CTSS.chr(CTSScoordinatesGR(ce))

.CTSS.chr <-
  setClass( "CTSS.chr"
          , contains = "CTSS"
          , validity =
  function(object) 
    if (length(unique(seqnames(object))) != 1)
      return("Mutiple sequnames found: CTSS.chr objects should be only on a single chromosome.")
)
