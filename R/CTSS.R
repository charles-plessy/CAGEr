#' CAGE Transcription Start Sites
#' 
#' The `CTSS` class represents CAGE transcription start sites (CTSS) at
#' single-nucleotide resolution, using [`GenomicRanges::UnstitchedGPos`] as base
#' class.  It is used by _CAGEr_ for type safety.
#' 
#' The `genomeName` element of the `metadata` slot is used to store the
#' name of the _BSgenome_ package used when constructing the `CAGEr` object.
#' 
#' @author Charles Plessy
#' 
#' @import methods
#' @importFrom GenomicRanges GRanges GPos from_GPos_to_GRanges
#' @export

setClass("CTSS", contains = "UnstitchedGPos")

#' @rdname CTSS-class
#' @param object See  [`methods::show`]

setMethod("show", "CTSS", function(object) {
  callNextMethod()
  bsgenomeName  <- ifelse(is.null(object@metadata$genomeName),
                          "unspecified",
                          object@metadata$genomeName)
  cat("  BSgenome name:", bsgenomeName , "\n")
})


#' @rdname CTSS-class
#' @param .Object See  [`methods::new`]
#' @examples 
#' # Convert an UnstitchedGPos object using the new() constructor.
#' gp <- GPos("chr1:2:-", stitch = FALSE)
#' ctss <- new("CTSS", gp, bsgenomeName = "BSgenome.Drerio.UCSC.danRer7")
#' genomeName(ctss)

setMethod("initialize", "CTSS", function(.Object, ..., bsgenomeName = NULL) {
  .Object <- callNextMethod(.Object, ...)
  if (! is.null(bsgenomeName))
    metadata(.Object)$genomeName <- bsgenomeName
  .Object
})


#' CTSS
#' 
#' The `CTSS` constructor takes the same arguments as [`GenomicRanges::GPos`],
#' plus `bsgenomeName`, and minus `stitch`, which is hardcoded to `FALSE`.
#' 
#' @param seqnames,pos,strand,seqinfo,seqlengths,... See the documentation
#'        of [`GenomicRanges::GPos`] for further details.
#' @param bsgenomeName String containing the name of a _BSgenome_ package.
#' 
#' @rdname CTSS-class
#' @examples 
#' 
#' # Create a new object using the CTSS() constructor.
#' CTSS("chr1", 2, "-", bsgenomeName = "BSgenome.Drerio.UCSC.danRer7")
#' @export

CTSS <- function (seqnames = NULL, pos = NULL, strand = NULL, ...,
                  seqinfo = NULL, seqlengths = NULL,
                  bsgenomeName = NULL) {
  gpos <- GPos(seqnames = seqnames, pos = pos, strand = strand,
               ..., seqinfo = seqinfo, seqlengths = seqlengths,
               stitch=FALSE)
  new("CTSS", gpos, bsgenomeName = bsgenomeName)
}

#' @rdname CTSS-class
#' @examples 
#' 
#' # Coerce CTSS to GRanges
#' as(ctss, "GRanges")
#' @export
# See https://stat.ethz.ch/pipermail/bioc-devel/2019-September/015524.html

setMethod("coerce", c("CTSS", "GRanges"), from_GPos_to_GRanges)

#' @rdname CTSS-class
#' 
#' @param from,to,strict See  [`methods::coerce`].
#' 
#' @details Coercion from `GRanges` to `CTSS` loses information, but it seems
#' to be fine, since other coercions like `as(1.2, "integer")` do the same.
#' 
#' @examples
#' 
#' # Coerce a GRanges object to CTSS using the as() method.
#' gr <- GRanges("chr1:1-10:-")
#' gr$seq <- "AAAAAAAAAA"
#' seqlengths(gr) <- 100
#' genome(gr) <- "foo"
#' as(gr, "CTSS")
#' identical(seqinfo(gr), seqinfo(as(gr, "CTSS")))
#' as(as(gr, "CTSS"), "CTSS") # Make sure it works twice in a row
#' @export

setMethod("coerce", c("GRanges", "CTSS"),
          function(from, to = "CTSS", strict = TRUE) {
  gr <- promoters(from, 0, 1)
  gp <- GPos(gr, stitch = FALSE)
  mcols(gp) <- mcols(gr)
  as(gp, "CTSS")
})


#' @description The `CTSS.chr` class represents a `CTSS` object that is 
#' guaranteed to be only on a single chromosome.  It is used internally by
#' _CAGEr_ for type-safe polymorphic dispatch.
#' 
#' @rdname CTSS-class
#' 
#' @import methods
#' @importFrom GenomicRanges GPos
#' 
#' @examples
#' 
#' # (internal use) Transform CTSS to CTSS.chr object
#' ctss.chr <- as(CTSScoordinatesGR(exampleCAGEexp), "CTSS.chr")

setClass( "CTSS.chr"
        , contains = "CTSS"
        , validity =
  function(object) 
    if (length(seqlevelsInUse(object)) > 1)
      return("Mutiple sequnames found: CTSS.chr objects should be only on a single chromosome.")
)


#' ConsensusClusters
#' 
#' The \code{ConsensusClusters} class represents consensus clusters.
#' It is used internally by CAGEr for type safety.
#' 
#' Consensus clusters must not overlap, so that a single TSS in the
#' genome can only be attributed to a single cluster.
#' @aliases ConsensusClusters

.ConsensusClusters <-
  setClass( "ConsensusClusters"
          , contains = "GRanges"
          , validity = function(object)
            if (length(reduce(object, min.gapwidth=0L)) < length(object))
              return("Consensus clusters must not overlap with each other."))

#' TagClusters
#'  
#' @details The \code{TagClusters} class represents tag clusters.
#' It is used internally by CAGEr for type safety.
#' @aliases TagClusters

.TagClusters <- setClass( "TagClusters", contains = "GRanges")
