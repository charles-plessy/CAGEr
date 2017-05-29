######################################################################################
# Funtions for normalizing CAGE tag count to a referent power-law distribution
# Reference: Balwierz, P. J., Carninci, P., Daub, C. O., Kawai, J., Hayashizaki, 
# Y., Van Belle, W., Beisel, C., et al. (2009). Methods for analyzing deep sequencing 
# expression data: constructing the human and mouse promoterome with deepCAGE data. 
# Genome Biology, 10(7), R79.
#

#' normalizeTagCount
#' 
#' Normalise CTSS tag counts in a CAGEr object.
#' 
#' @param object A CAGEr object of the CAGEset or CAGEexp class.
#' 
#' @return A CAGEr object of the same class, containing norm datalised data
#' with the \code{CTSSnormalizedTpm()} functions.
#' 
#' @family CAGEr normalised data functions
#' 
#' @seealso CTSSnormalizedTpm, CTSSnormalizedTpmDf

setGeneric(
name="normalizeTagCount",
def=function(object, method = "powerLaw", fitInRange = c(10, 1000), alpha = 1.25, T = 10^6){
	standardGeneric("normalizeTagCount")
}
)

#' .normalizeTagCount_switcher
#' 
#' Common code to normalizeTagCount for CAGEset and CAGEexp objects
#' Do not reuse elsewhere.
#' 
#' @param method The method.
#' @param object A CAGEset or CAGEexp object.
#' 
#' @return A data.frame for CAGEset objects, or a DataFrame for CAGEexp objects.

.normalizeTagCount_switcher <- function(method, object, fitInRange, alpha, T) {
  message("\nNormalizing tag count...")
  switch( method
        , powerLaw  = CAGEr:::.powerLaw(CTSStagCountTable(object), fitInRange, alpha, T)
        , simpleTpm = CAGEr:::.simpleTpm(CTSStagCountTable(object))
        , none      = CTSStagCountTable(object)
        , stop('"method" must be one of ("powerLaw", "simpleTpm", "none")'))
}

# For the CAGEset class, normalizeTagCount populates the normalizedTpmMatrix slot.
 
setMethod("normalizeTagCount",
signature(object = "CAGEset"),
function (object, method, fitInRange, alpha, T){
	objName <- deparse(substitute(object))
	object@normalizedTpmMatrix <- .normalizeTagCount_switcher(method, object, fitInRange, alpha, T)
	assign(objName, object, envir = parent.frame())
	invisible(1)
})

# For the CAGEexp class, normalizeTagCount populates the normalized slot or the tagCountMatrix
# experiment.

setMethod("normalizeTagCount",
signature(object = "CAGEexp"),
function (object, method, fitInRange, alpha, T){
	objName <- deparse(substitute(object))
	se <- CTSStagCountSE(object)
	assays(se)$normalizedTpmMatrix <- .normalizeTagCount_switcher(method, object, fitInRange, alpha, T)
  CTSStagCountSE(object) <- se
	assign(objName, object, envir = parent.frame())
	invisible(1)
})
