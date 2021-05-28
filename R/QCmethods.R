#' @include CAGEr.R

#' @name summariseChrExpr
#' @title Expression levels by chromosomes
#' 
#' @description Counts the number of molecules detected per chromosome, normalises
#' by library size and stores the raw and normalised results in the
#' \code{\link{CAGEr}} object.
#' 
#' @param object A `CAGEexp` object objects are not supported).
#' 
#' @return Modifies the \code{CAGEexp} by adding a \dQuote{seqNameTotals} experiment
#' containing matrices where rows represent chromosomes and columns represent samples.
#' 
#' @author Charles Plessy
#' 
#' @family CAGEr object modifiers
#' @seealso seqNameTotals
#' 
#' @examples 
#' summariseChrExpr(exampleCAGEexp)
#' 
#' @export

setGeneric("summariseChrExpr", function(object) standardGeneric("summariseChrExpr"))

#' @rdname summariseChrExpr

setMethod("summariseChrExpr", "CAGEexp", function(object) {
  objname <- deparse(substitute(object))
  seqNameTotals <- lapply( sampleList(object)
                         , function(n) tapply( CTSStagCountDF(object)[,n]
                                             , seqnames(CTSScoordinatesGR(object))
                                             , sum))
  seqNameTotals <- as.matrix(as.data.frame(seqNameTotals))
  seqNameTotals[is.na(seqNameTotals)] <- 0
  seqNameTotals.norm <- prop.table(seqNameTotals, margin = 2) * 100
  seqNameTotalsSE(object) <-
    SummarizedExperiment( assays  = SimpleList( counts = seqNameTotals
                                              , norm   = seqNameTotals.norm)
                        , rowData = DataFrame(as.data.frame(seqinfo(getRefGenome(genomeName(object))))))
  validObject(object)
  object
})
