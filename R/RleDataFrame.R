#' rowSums function for Rle DataFrames
#' 
#' Drop-in replacement for the `rowSums` function, which does not work natively
#' on [`S4Vectors::DataFrame`] objects containing [`S4Vectors::Rle`]-encoded
#' numerical values.
#' 
#' See the file `benchmarks/rowSums_on_Rle_DF.md` in the source Git repository
#' of _CAGEr_ for the alternatives that were considered.
#' 
#' @param x A `DataFrame` containing only numerical `Rle` columns.
#' 
#' @param na.rm logical. Should missing values (including `NaN`) be omitted from
#'        the calculations?
#'        
#' @return A `Rle`-encoded numerical vector of the same class as in the
#' `DataFrame`.
#' 
#' @author Charles Plessy
#' @family Rle DataFrames
#' 
#' @examples 
#' exampleCAGEexp |> CTSStagCountDF() |> CAGEr:::rowSums.RleDataFrame(na.rm = TRUE)

rowSums.RleDataFrame <- function (x, na.rm = FALSE) {
  if(isTRUE(na.rm)) x <- lapply(x, \(v) {v[is.na(v)] <- 0L; v})
  Rle(Reduce(`+`, lapply(x, decode), 0L))
}

#' rowsum function for Rle DataFrames
#' 
#' Drop-in replacement for the `rowsum` function, which does not work natively
#' on [`S4Vectors::DataFrame`] objects containing [`S4Vectors::Rle`]-encoded
#' numerical values.
#' 
#' See the file `benchmarks/rowsum_on_Rle_DF.md` in the source Git repository
#' of _CAGEr_ for the alternatives that were considered.
#' 
#' @param x A `DataFrame` containing only numerical `Rle` columns.
#' 
#' @param group a vector or factor giving the grouping, with one element per row
#'        of `x`.  Missing values will be treated as another group and a warning
#'        will be given.
#' 
#' @param reorder If `TRUE`, then the result will be in order of
#'        `sort(unique(group))`, if `FALSE`, it will be in the order that groups
#'        were encountered.
#' 
#' @param na.rm Logical (`TRUE` or `FALSE`). Should `NA` (including `NaN`)
#'        values be discarded?
#' 
#' @param ... Other arguments to be passed to or from methods.
#' 
#' @author Charles Plessy
#' @family Rle DataFrames
#' 
#' @examples 
#' exampleCAGEexp |> CTSStagCountDF() |>
#'   CAGEr:::rowsum.RleDataFrame(decode(CTSScoordinatesGR(exampleCAGEexp)$cluster), reorder = FALSE)

rowsum.RleDataFrame <- function(x, group, reorder = TRUE, na.rm = FALSE, ...)
  rowsum(
    as.data.frame(lapply(x, decode)),
    group,
    reorder = reorder,
    na.rm = FALSE,
    ...) |> DataFrame()
