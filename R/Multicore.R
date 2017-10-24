#' @name CAGEr_Multicore
#' 
#' @title Multicore support in CAGEr
#' 
#' @description CAGEr is in the transition towards using the BiocParallel for
#' multicore parallelisation.  On Windows platforms, the multicore support
#' is disabled transparently, that is, attempts to use multiple cores are
#' silently ignored.
#' 
#' @param useMulticore TRUE or FALSE
#' @param nrCores number of cores to use (leave \code{NULL} to let BiocParallel
#'        choose).
#' 
#' @return Returns either a \code{MulticoreParam} object or a
#' \code{SerialParam} object.
#' 
#' @author Charles Plessy
#' 
#' @examples 
#' CAGEr:::CAGEr_Multicore()
#' CAGEr:::CAGEr_Multicore(TRUE,)
#' CAGEr:::CAGEr_Multicore(TRUE, 3)
#' CAGEr:::CAGEr_Multicore(FALSE, 3)
#' 
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel MulticoreParam
#' @importFrom BiocParallel multicoreWorkers
#' @importFrom BiocParallel SerialParam

CAGEr_Multicore <- function (useMulticore = FALSE, nrCores = NULL) {
  if (useMulticore) {
    if (is.null(nrCores)) nrCores <- multicoreWorkers()
    MulticoreParam(workers = nrCores)
  }	else {
    SerialParam()
  }
}
