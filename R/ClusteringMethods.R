#' @include CTSS.R Multicore.R

#' @rdname byCtss
#' 
#' @title Apply functions to identical CTSSes.
#' 
#' @param ctssDT A [`data.table::data.table`] representing CTSSes.
#' @param colName The name of the column on which to apply the function.
#' @param fun The function to apply.
#' 
#' @description `.byCTSS` is a private function using  `data.table` objects
#' to preform grouping operations at a high performance.  These functions use
#' _non-standard evaluation_ in a context that raises warnings in `R CMD check`.
#' By separating these functions from the rest of the code, I hope to make the workarounds
#' easier to manage.
#' 
#' @examples
#' ctssDT <- data.table::data.table(
#'   chr       = c("chr1", "chr1", "chr1", "chr2"),
#'   pos       = c(1     , 1     , 2     , 1     ),
#'   strand    = c("+"   , "+"   , "-"   , "-"   ),
#'   tag_count = c(1     , 1     , 1     , 1     ))
#' ctssDT
#' CAGEr:::.byCtss(ctssDT, "tag_count", sum)

setGeneric( ".byCtss"
          , function (ctssDT, colName, fun) standardGeneric(".byCtss"))

#' @rdname byCtss

setMethod(".byCtss", "data.table", function (ctssDT, colName, fun) {
  if (! all(c("chr", "pos", "strand") %in% colnames(ctssDT))) stop("These are not CTSSes.")
  chr <- pos <- strand <- .SD <- NULL
  ctssDT[ , fun(.SD[[1]])
        , by = list(chr, pos, strand)
        , .SDcols = colName]
})
