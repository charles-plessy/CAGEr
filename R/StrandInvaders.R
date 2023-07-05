#' @rdname strandInvaders
#' 
#' @name Strand invaders
#' 
#' @title Detect and remove strand invasion artefacts
#' 
#' @description `findStrandInvaders` detects strand invasion artefacts in the
#' CTSS data.  `removeStrandInvaders` removes them.
#' 
#' _Strand invaders_ are artefacts produced by _template switching_ reactions
#' used in methods such as _nanoCAGE_ and its derivatives (_C1 CAGE_, ...).
#' They are described in details in Tang _et al._, 2013.  Briefly, these
#' artefacts create CAGE-like signal downstream of genome sequences highly
#' similar to the tail of template-switching oligonucleotides, which is
#' `TATAGGG` in recent (2017) nanoCAGE protocols.  Since these artefacts
#' represent truncated cDNAs, they do not indicate promoter regions.  It is
#' therefore advisable to remove these artefacts.  Moreover, when a sample
#' barcode is near the linker sequence (which is not the case in recent
#' nanoCAGE protocols), the strand-invasion artefacts can produce
#' _sample-specific biases_, which can be confounded with biological effects
#' depending on how the barcode sequences were chosen.  A `barcode` parameter
#' is provided to incorporate this information.
#' 
#' @param object A [`CAGEexp`] object object containing CTSS data and the name
#'        of a reference genome.
#'        
#' @param distance The maximal edit distance between the genome and linker
#'        sequences.  Regardless this parameter, only a single mismatch is
#'        allowed in the last three bases of the linker.
#' 
#' @param barcode A vector of sample barcode sequences, or the name of a column
#'        metadata of the `CAGEexp` object containing this information.
#'        (_Not implemented yet_)
#' 
#' @param linker The sequence of the tail of the template-switching
#'        oligonucleotide, that will be matched with the genome sequence
#'        (defaults to `TATAGGG`).
#' 
#' @return `findStrandInvaders` returns a logical-[Rle] vector indicating the
#' position of the strand invaders in the input ranges.
#' 
#' @return With [CTSS] objects as input `removeStrandInvaders` returns the
#' object after removing the CTSS positions identified as strand invaders.
#' In the case of `CAGEexp` objects, a modified object is returned.  Its sample
#' metadata is also updated by creating a new `strandInvaders` column that
#' indicates the number of molecule counts removed.  This value is subtracted
#' from the `counts` colum so that the total number of tags is still equal to
#' `librarySizes`.
#' 
# @family CAGEr object modifiers
# @family CAGEr TSS functions
#' 
#' @references Tang _et al._, \dQuote{Suppression of artifacts and barcode bias in
#' high-throughput transcriptome analyses utilizing template switching.}
#' _Nucleic Acids Res._ **2013** Feb 1;41(3):e44.
#' PubMed ID: [23180801](https://pubmed.gov/23180801),
#' DOI: [10.1093/nar/gks112](https://doi.org/10.1093/nar/gks1128)
#' 
#' @examples 
#' # Note that these examples do not do much on the example data since it was
#' # not constructed using a protocol based using the template-switching method.
#' 
#' findStrandInvaders(exampleCAGEexp)
#' removeStrandInvaders(exampleCAGEexp)
#' 

NULL

#' @rdname strandInvaders
#' @export

setGeneric( "findStrandInvaders"
          , function( object, distance = 1, barcode = NULL, linker = "TATAGGG")
              standardGeneric("findStrandInvaders"))

#' @rdname strandInvaders
#' @export

setGeneric( "removeStrandInvaders"
          , function( object, distance = 1, barcode = NULL, linker = "TATAGGG")
              standardGeneric("removeStrandInvaders"))

#' @rdname strandInvaders
#' @export

setMethod( "findStrandInvaders", "CAGEexp"
         , function (object, distance, barcode, linker)
  findStrandInvaders(CTSScoordinatesGR(object), distance, barcode, linker))

#' @rdname strandInvaders
#' @export

setMethod( "removeStrandInvaders", "CAGEexp"
         , function (object, distance, barcode, linker) {
  if ( ! is.null(barcode)) stop("Sorry, the barcode is not implemented yet.\nPlease ask by opening an issue on GitHub.")
  strandInvaders <- findStrandInvaders(object, distance, barcode, linker)
  se <- CTSStagCountSE(object)
  CTSStagCountSE(object) <- se[!decode(strandInvaders),]
  newcounts <- sapply(CTSStagCountDF(object), sum)
  object$strandInvaders <- object$librarySizes - newcounts
  object$librarySizes <- newcounts
  object
})

#' @rdname strandInvaders
#' @export

setMethod( "findStrandInvaders", "CTSS"
         , function (object, distance, barcode, linker) {
  .flagByUpstreamSequences(object, target = linker, distance = distance, strandInvaders = TRUE)
})

#' @rdname strandInvaders
#' @export

setMethod( "removeStrandInvaders", "CTSS"
         , function (object, distance, barcode, linker)
  object[!decode(findStrandInvaders(object))])


#' Filter by upstream sequences
#' 
#' Looks up the bases directly upstream provided _genomic ranges_ and searches
#' for a gapless match with a _target_ seqence within a given edit _distance_.
#' 
#' If the provided `object` represents _tag clusters_ or _consensus clusters_,
#' the search will be done upstream its _dominant peak_.  Convert the object
#' to the `GRanges` class if this is not the behaviour you want.
#' 
#' @param object A [`CTSS`], a [`TagClusters`], [`ConsensusClusters`] or a
#'        [`GenomicRanges::GRanges`] object from which a _BSgenome_ object can
#'        be reached.
#' 
#' @param target A target sequence.
#' 
#' @param distance The maximal edit distance between the genome and the target
#'        sequence (default: 0).
#' 
#' @returns A `logical-RLe` vector indicating if ranges matched the target.
#' 
#' @importFrom BSgenome getBSgenome
#' @importFrom BSgenome getSeq
#' @importFrom stringdist stringdist
#' @importFrom stringi stri_sub
#' 
#' @author Charles Plessy
#' 
#' @family CAGEr filter functions
#' 
#' @export

setGeneric( "flagByUpstreamSequences"
          , function (object, target, distance = 0)
            standardGeneric("flagByUpstreamSequences"))

#' @rdname flagByUpstreamSequences
#' @export

setMethod( "flagByUpstreamSequences", "CTSS"
         , function (object, target, distance = 0) {
  flagByUpstreamSequences(GRanges(object), target = target, distance = distance)
})

#' @rdname flagByUpstreamSequences
#' @export

setMethod( "flagByUpstreamSequences", "TagClusters"
         , function (object, target, distance = 0) {
  flagByUpstreamSequences(object$dominant_ctss, target = target, distance = distance)
})

#' @rdname flagByUpstreamSequences
#' @export

setMethod( "flagByUpstreamSequences", "ConsensusClusters"
         , function (object, target, distance = 0) {
  flagByUpstreamSequences(object$dominant_ctss, target = target, distance = distance)
})

#' @rdname flagByUpstreamSequences
#' @export

setMethod( "flagByUpstreamSequences", "GRanges"
           , function (object, target, distance = 0) {
 .flagByUpstreamSequences(GRanges(object), target = target, distance = distance)
})

#' @noRd
#' @examples 
#' # Note that we do not expect TATAGGG artefacts in this example dataset
#' CAGEr:::.flagByUpstreamSequences(CTSScoordinatesGR(exampleCAGEexp), "TATAGGG", 2)
#' CAGEr:::.flagByUpstreamSequences(CTSScoordinatesGR(exampleCAGEexp), "TATAGGG", 2, strandInvaders = TRUE)
           
.flagByUpstreamSequences <- function (object, target, distance = 0, strandInvaders = FALSE) {
  bsgenome <- object |> genome() |> unique() |> getBSgenome()
  updreamRanges <- promoters(object, nchar(target), 0) |> trim() |> suppressWarnings()
  upstreamSeq <- getSeq(bsgenome, updreamRanges)
  matching <- Rle(stringdist(upstreamSeq, target) <= distance)
  # Only remove if at least 2 of the last 3 nucleotides are guanosines,
  # as in https://github.com/davetang/23180801/blob/master/find_strand_invasion-20130307.pl
  if (distance > 1 & isTRUE(strandInvaders)) {
   matching <-
     matching &
      stringdist( stri_sub(upstreamSeq, nchar(target) - 2, nchar(target))
                , stri_sub(target, nchar(target) - 2, nchar(target)))       < 2
  }
  matching
}
