######################################################################################
# Funtions for normalizing CAGE tag count to a referent power-law distribution
# Reference: Balwierz, P. J., Carninci, P., Daub, C. O., Kawai, J., Hayashizaki, 
# Y., Van Belle, W., Beisel, C., et al. (2009). Methods for analyzing deep sequencing 
# expression data: constructing the human and mouse promoterome with deepCAGE data. 
# Genome Biology, 10(7), R79.
#

#' @name normalizeTagCount
#' 
#' @title Normalizing raw CAGE tag count
#' 
#' @description Normalizes raw CAGE tag count per CTSS in all experiments to a same
#' referent distribution.  A simple tag per million normalization or normalization to a
#' referent power-law distribution (Balwierz \emph{et al}., Genome Biology 2009) can be
#' specified.
#' 
#' @param object A \code{\link{CAGEset}} object
#' 
#' @param method Method to be used for normalization.  Can be either \code{"simpleTpm"}
#'   to convert tag counts to tags per million or \code{"powerLaw"} to normalize to a
#'   referent power-law distribution, or \code{"none"} to keep using the raw tag counts
#'   in downstream analyses.
#' 
#' @param fitInRange An integer vector with two values specifying a range of tag count
#'   values to be used for fitting a power-law distribution to reverse cumulatives.
#'   Used only when \code{method = "powerLaw"}, otherwise ignored.  See Details.
#' 
#' @param alpha	\code{-1 * alpha} will be the slope of the referent power-law distribution
#'   in the log-log representation.  Used only when \code{method = "powerLaw"}, otherwise
#'   ignored.  See Details.
#'   
#' @param T Total number of CAGE tags in the referent power-law distribution.  Setting
#'   \code{T = 10^6} results in normalized values that correspond to tags per million in
#'   the referent distribution.  Used only when \code{method = "powerLaw"}, otherwise
#'   ignored.  See Details.
#' 
#' @details It has been shown that many CAGE datasets follow a power-law distribution
#' (Balwierz \emph{et al}., Genome Biology 2009).  Plotting the number of CAGE tags
#' (X-axis) against the number of TSSs that are supported by >= of that number of tags
#' (Y-axis) results in a distribution that can be approximated by a power-law.  On a
#' log-log scale this theoretical referent distribution can be described by a
#' monotonically decreasing linear function \code{y = -1 * alpha * x + beta}, which is
#' fully determined by the slope \code{alpha} and total number of tags \code{T} (which
#' together with \code{alpha} determines the value of \code{beta}).  Thus, by specifying
#' parameters \code{alpha} and \code{T} a desired referent power-law distribution can be
#' selected.  However, real CAGE datasets deviate from the power-law in the areas of very
#' low and very high number of tags, so it is advisable to discard these areas before
#' fitting a power-law distribution.  \code{fitInRange} parameter allows to specify a
#' range of values (lower and upper limit of the number of CAGE tags) that will be used to
#' fit a power-law.  Plotting reverse cumulatives using \code{\link{plotReverseCumulatives}}
#' function can help in choosing the best range of values.  After fitting a power-law
#' distribution to each CAGE dataset individually, all datasets are normalized to a
#' referent distribution specified by \code{alpha} and \code{T}. When \code{T = 10^6},
#' normalized values are expressed as tags per million (tpm).
#' 
#' @return The slot \code{normalizedTpmMatrix} of the provided \code{\link{CAGEset}}
#' object will be occupied by normalized CAGE signal values per CTSS across all experiments,
#' or with the raw tag counts (in case \code{method = "none"}).
#' 
#' @references
#' 
#' Balwierz \emph{et al}. (2009) Methods for analyzing deep sequencing expression data:
#' constructing the human and mouse promoterome with deepCAGE data, \emph{Genome Biology}
#' \bold{10}(7):R79. 
#' 
#' @author Vanja Haberle
#' 
#' @seealso \code{\link{plotReverseCumulatives}}, \code{\link{CTSSnormalizedTpm}}
#' 
#' @family CAGEr object modifiers
#' @family CAGEr normalised data functions
#' 
#' @examples 
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' normalizeTagCount(exampleCAGEset, method = "powerLaw")
#' 
#' @export

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
#' 
#' @noRd

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
