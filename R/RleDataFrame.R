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
#' 
#' @examples 
#' exampleCAGEexp |> CTSStagCountDF() |> rowSums.RleDataFrame(na.rm = TRUE)

rowSums.RleDataFrame <- function (x, na.rm = FALSE) {
  if(isTRUE(na.rm)) x <- lapply(x, \(v) {v[is.na(v)] <- 0L; v})
  Rle(Reduce(`+`, lapply(x, decode), 0L))
}
