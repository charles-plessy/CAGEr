#' TSS logo
#' 
#' Plot the sequence logo of the region flanking the TSS.  When this function
#' is given _tag clusters_ or _consensus clusters_, it uses the _dominant peak_
#' as the transcription start site.
#' 
#' This function will only work if the [`CAGEexp`] object was built with a
#' [`BSgenome`] package, as it needs to extract genomic sequences.
#' 
#' @param x A [`CTSS`], a [`TagClusters`] or a [`ConsensusClusters`] object.
#' @param upstream Number of bases to plot upstream the TSS.
#' @param downstream Number of bases to plot downstream the TSS, including the
#'        TSS itself.
#' 
#' @returns A [`ggplot2::ggplot`] object showing the sequence logo.  The
#' coordinates displayed are negative for upstream sequences and positive
#' downstream.  The position of the TSS is set to 1.
#' 
#' @family CAGEr plot functions
#' @family CAGEr TSS functions
#' 
#' @author Charles Plessy
#' 
#' @examples
#' TSSlogo(exampleCAGEexp|>consensusClustersGR(), 20, 10)
#' 
#' @importFrom Biostrings consensusMatrix
#' 
#' @export

setGeneric("TSSlogo",
           function (x, upstream = 10, downstream = 10)
             standardGeneric("TSSlogo"))

#' @rdname TSSlogo

setMethod("TSSlogo", signature("CAGEexp"),
          function (x, upstream = 10, downstream = 10) {
  TSSlogo(CTSScoordinatesGR(object), upstream = upstream, downstream = downstream)
})

#' @rdname TSSlogo

setMethod("TSSlogo", signature("TagClusters"),
          function (x, upstream = 10, downstream = 10) {
  TSSlogo(object$dominant_ctss |> as("CTSS"), upstream = upstream, downstream = downstream)
})

setMethod("TSSlogo", signature("ConsensusClusters"),
          function (x, upstream = 10, downstream = 10) {
  TSSlogo(object$dominant_ctss |> as("CTSS"), upstream = upstream, downstream = downstream)
})

#' @rdname TSSlogo

setMethod("TSSlogo", signature("CTSS"),
          function (x, upstream = 10, downstream = 10) {
  .TSSlogo(object, upstream = upstream, downstream = downstream)
})

.TSSlogo <- function(x, upstream=10, downstream=10) {
  if (! requireNamespace("ggseqlogo"))
    stop("This function requires the ", dQuote("ggseqlogo"), " package; please install it.")
  # Extract sequences
  upstreamRanges <-
    promoters(x = x, upstream = upstream, downstream = downstream) |>
    suppressWarnings() # This warns about off-genome coordinates, but we will fix this below.
  
  # Discard ranges that would be trimmed 
  upstreamRanges_trimmed <- upstreamRanges |> trim()
  upstreamRanges <- upstreamRanges[width(upstreamRanges) == width(upstreamRanges_trimmed)]
  
  # Extract sequences
  bsgenome <- getBSgenome(unique(genome(x)))
  upstreamSeq <- getSeq(bsgenome, upstreamRanges)
  
  # Plot sequence logo
  letter_counts <- consensusMatrix(upstreamSeq)
  probs <- prop.table(letter_counts[1:4,], 2)
  gg <- ggseqlogo::ggseqlogo(probs)
  # Circumvent "Scale for x is already present." warning.
  gg$scales$scales[[2]]$breaks <- seq(    1       , upstream + downstream, by = 1)
  gg$scales$scales[[2]]$labels <- seq(1 - upstream,            downstream, by = 1)
  gg$scales$scales[[1]]$breaks <- seq(    1       , upstream + downstream, by = 1)
  gg$scales$scales[[1]]$labels <- seq(1 - upstream,            downstream, by = 1)
  gg
}