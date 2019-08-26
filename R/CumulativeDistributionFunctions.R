#' @name coverage-functions
#' @title Private functions behind `cumulativeCTSSdistribution`
#' @param clusters GRanges as per [tagClustersGR()].
#' @param useMulticore,nrCores See clusterCTSS.
#' @examples 
#' library(GenomicRanges)
#' library(IRanges)
#' ctss <- CAGEr:::.CTSS(
#'           GPos( stitch = FALSE
#'               , GRanges( seqnames=Rle("chr1")
#'                        , IRanges(c(1,3,4,12,14,25,28,31,35), w=1)
#'                        , strand = "+")))
#' score(ctss) <- 1
#' ctss.chr <- CAGEr:::.CTSS.chr(ctss)
#' clusters <- GRanges( seqnames = Rle("chr1")
#'                    , ranges = IRanges(c(1,12,25,31,32), c(4,14,28,31,33))
#'                    , strand = "+")
#' chrom <- "chr1"
#' str <- "+"
NULL

#' @name .getCAGEsignalCoverage
#' @rdname coverage-functions
#' 
#' @details `.getCAGEsignalCoverage`` does...
#' Note that strand is not taken into account.
#' 
#' @param ctss.chr A `CTSS.chr`` object (guaranteed to have only one chromosome).
#' 
#' @importFrom IRanges Views
#' @importFrom IRanges viewApply
#' 
#' @examples
#' CAGEr:::.getCAGEsignalCoverage(ctss.chr, clusters)

setGeneric(".getCAGEsignalCoverage", function(ctss.chr, clusters) standardGeneric(".getCAGEsignalCoverage"))

setMethod(".getCAGEsignalCoverage", c("CTSS.chr", "GRanges"), function(ctss.chr, clusters) {
  cov <- Rle(rep(0, max(end(clusters), end(ctss.chr))))
	cov[start(ctss.chr)] <- score(ctss.chr)
  cluster.cumsums <- Views(cov, start = start(clusters), end = end(clusters))
	viewApply(cluster.cumsums, cumsum)
})


#' @name .getCumsumChr2
#' @rdname coverage-functions
#' 
#' @details `.getCumsumChr2`
#' 
#' @param chrom a chromosome name
#' @param str a strand name
#' 
#' @examples 
#' CAGEr:::.getCumsumChr2(clusters, ctss, chrom, str)

setGeneric(".getCumsumChr2", function(clusters, ctss, chrom, str) standardGeneric(".getCumsumChr2"))

setMethod(".getCumsumChr2", c("GRanges", "CTSS"), function(clusters, ctss, chrom, str) {
  clusters.chr <- clusters[seqnames(clusters) == chrom & strand(clusters) == str]
  ctss.chr <- ctss[seqnames(ctss) == chrom & strand(ctss) == str]
  if (length(clusters.chr) > 0 & length(ctss) > 0) {
    .getCAGEsignalCoverage(ctss.chr = .CTSS.chr(ctss.chr), clusters = clusters.chr)
  } else {
    return()
  }
})


#' @name .getCumsum
#' @rdname coverage-functions
#' 
#' @description `.getCumsum` calculates cumulative sums of tpm along the clusters.
#' 
#' @param ctss GRanges as per `CTSScoordinatesGR`, with the score of one sample.
#' 
#' @return `.getCumsum` returns a list of `Rle` vectors (IRanges package) containing cumulative
#' sum for each cluster (length of list is equal to number of clusters and names of the list
#' components corespond to the name of the corresponding cluster) v.
#'
#' @examples 
#' ctss      <- CTSSnormalizedTpmGR(exampleCAGEset, "sample1")
#' ctss      <- ctss[ctss$filteredCTSSidx]
#' clusters  <- tagClustersGR(exampleCAGEset, "sample1")
#' clusters.cumsum <- RleList(CAGEr:::.getCumsum(ctss, clusters))
#' identical( lapply(exampleCAGEset@CTSScumulativesTagClusters[[1]],decode)
#'          , lapply(clusters.cumsum, decode))
#' # Not identical if not decoded because Rle method is attached to S4Vectors in one case
#' # and to IRanges in the other case.
#' decode(clusters.cumsum[[1]])
#' ctss[queryHits(findOverlaps(ctss, clusters[1]))]
#' clusters[1]
#' 
#' ctss      <- CTSSnormalizedTpmGR(exampleCAGEexp, "Zf.30p.dome")
#' ctss      <- ctss[ctss$filteredCTSSidx]
#' clusters  <- tagClustersGR(exampleCAGEexp, "Zf.30p.dome")
#' clusters.cumsum <- CAGEr:::.getCumsum(ctss, head(clusters))
#' decode(clusters.cumsum[[3]])
#' ctss[queryHits(findOverlaps(ctss, clusters[3]))]
#' clusters[3]

setGeneric(".getCumsum", function(ctss, clusters, useMulticore = FALSE, nrCores = NULL) standardGeneric(".getCumsum"))

setMethod(".getCumsum", c("CTSS", "GRanges"), function(ctss, clusters, useMulticore , nrCores) {
  getCumSumChrStrand <- function(chrom) {
    plus.cumsum  <- .getCumsumChr2(clusters = clusters, ctss = ctss, chrom = chrom, str = "+")
    minus.cumsum <- .getCumsumChr2(clusters = clusters, ctss = ctss, chrom = chrom, str = "-")
    c(plus.cumsum, minus.cumsum)
  }
  clusters.cumsum <- unlist(bplapply( seqlevels(clusters), getCumSumChrStrand
                                    , BPPARAM = CAGEr_Multicore(useMulticore, nrCores)))
	names(clusters.cumsum) <- names(clusters)
	clusters.cumsum
})
