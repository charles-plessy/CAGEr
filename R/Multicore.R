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
#' @examples 
#' useMulticore <- CAGEr:::.checkMulticore(useMulticore)

.checkMulticore <- function (useMulticore) {
	pt <- .Platform$OS.type
	if(useMulticore == TRUE){
		if(pt == "unix"){
			if("parallel" %in% rownames(installed.packages()) == FALSE){
				stop("Cannot use multicore because package 'parallel' is not installed!")
			}
		}else{
			useMulticore = FALSE
			warning("Multicore is not supported on non-Unix platforms! Setting `useMulticore`` to FALSE")
		}
	}
	useMulticore
}
