#' @name CAGEr_Multicore
#' 
#' @title Multicore support in CAGEr
#' 
#' @description CAGEr is in the transition towards using the BiocParallel for
#' multicore parallelisation.  On Windows platforms, the multicore support
#' is disabled transparently, that is, attempts to use multiple cores are
#' silently ignored.
#' 
#' @author Charles Plessy
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

#' .checkMulticore
#' 
#' Check if multicore exectuion will be possible.
#' 
#' Non-exported private helper function.
#' 
#' @param useMulticore TRUE or FALSE
#' 
#' @return TRUE or FALSE
#' 
#' @family CAGEr multicore-enabled functions
#' 
#' @noRd
#' 
#' @importFrom utils installed.packages
#' @importFrom parallel mclapply
#' 
#' @examples 
#' useMulticore <- CAGEr:::.checkMulticore(useMulticore)

.checkMulticore <- function (useMulticore) {
	pt <- .Platform$OS.type
	if(useMulticore == TRUE){
		if(pt == "unix"){
			if("parallel" %in% rownames(installed.packages()) == FALSE){
				stop("Cannot use multicore because package 'parallel' is not installed!")
			}
		  requireNamespace("parallel")
		}else{
			useMulticore = FALSE
			warning("Multicore is not supported on non-Unix platforms! Setting `useMulticore`` to FALSE")
		}
	}
	useMulticore
}

.getNrCores <- function(nrCores)
  ifelse(is.null(nrCores), detectCores(), nrCores)