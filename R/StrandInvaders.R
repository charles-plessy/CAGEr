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
#' position of the strand invaders in the input ranges.  With [CTSS] objects as
#' input `removeStrandInvaders` returns the object after removing the CTSS
#' positions identified as strand invaders.  In the case of `CAGEexp` objects,
#' the input object is modified in place.  Its sample metadata is also
#' updated by creating a new `strandInvaders` column that indicates the number
#' of molecule counts removed.  This value is subtracted from the `counts` colum
#' so that the total number of tags is still equal to `librarySizes`.
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
#' "BSgenome.Drerio.UCSC.danRer7"
#' 
#' findStrandInvaders(exampleCAGEexp)
#' removeStrandInvaders(exampleCAGEexp)
#' 
#' @importFrom BSgenome getSeq
#' @importFrom stringdist stringdist
#' @importFrom stringi stri_sub
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
  upstreamSeq <- getSeq( getRefGenome(genomeName(object))
                       , suppressWarnings(trim(promoters(GRanges(object), nchar(linker), 0)))
                       , as.character = TRUE)
  strandInvaders <- Rle(stringdist(upstreamSeq, linker) <= distance)
  # Only remove if at least 2 of the last 3 nucleotides are guanosines,
  # as in https://github.com/davetang/23180801/blob/master/find_strand_invasion-20130307.pl
  if (distance > 1) {
  strandInvaders <-
    strandInvaders &
    stringdist( stri_sub(upstreamSeq, nchar(linker) - 2, nchar(linker))
              , stri_sub(linker, nchar(linker) - 2, nchar(linker))) < 2
  }
  strandInvaders
})

#' @rdname strandInvaders
#' @export

setMethod( "removeStrandInvaders", "CTSS"
         , function (object, distance, barcode, linker)
  object[!decode(findStrandInvaders(object))])
