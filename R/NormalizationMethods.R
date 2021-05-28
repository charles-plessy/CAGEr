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
#' referent power-law distribution (Balwierz _et al_., Genome Biology 2009) can be
#' specified.
#' 
#' @param object A [`CAGEexp`] object
#' 
#' @param method Method to be used for normalization.  Can be either `"simpleTpm"`
#'   to convert tag counts to tags per million or `"powerLaw"` to normalize to a
#'   referent power-law distribution, or `"none"` to keep using the raw tag counts
#'   in downstream analyses.
#' 
#' @param fitInRange An integer vector with two values specifying a range of tag count
#'   values to be used for fitting a power-law distribution to reverse cumulatives.
#'   Used only when `method = "powerLaw"`, otherwise ignored.  See Details.
#' 
#' @param alpha	\code{-1 * alpha} will be the slope of the referent power-law distribution
#'   in the log-log representation.  Used only when `method = "powerLaw"`, otherwise
#'   ignored.  See Details.
#'   
#' @param T Total number of CAGE tags in the referent power-law distribution.  Setting
#'   \code{T = 10^6} results in normalized values that correspond to tags per million in
#'   the referent distribution.  Used only when `method = "powerLaw"`, otherwise
#'   ignored.  See Details.
#' 
#' @details It has been shown that many CAGE datasets follow a power-law distribution
#' (Balwierz _et al_., Genome Biology 2009).  Plotting the number of CAGE tags
#' (X-axis) against the number of TSSs that are supported by >= of that number of tags
#' (Y-axis) results in a distribution that can be approximated by a power-law.  On a
#' log-log scale this theoretical referent distribution can be described by a
#' monotonically decreasing linear function `y = -1 * alpha * x + beta`, which is
#' fully determined by the slope `alpha` and total number of tags `T` (which
#' together with `alpha` determines the value of `beta`).  Thus, by specifying
#' parameters `alpha` and `T` a desired referent power-law distribution can be
#' selected.  However, real CAGE datasets deviate from the power-law in the areas of very
#' low and very high number of tags, so it is advisable to discard these areas before
#' fitting a power-law distribution. `fitInRange` parameter allows to specify a
#' range of values (lower and upper limit of the number of CAGE tags) that will be used to
#' fit a power-law.  Plotting reverse cumulatives using [`plotReverseCumulatives`]
#' function can help in choosing the best range of values.  After fitting a power-law
#' distribution to each CAGE dataset individually, all datasets are normalized to a
#' referent distribution specified by `alpha` and `T`. When `T = 10^6`,
#' normalized values are expressed as tags per million (tpm).
#' 
#' @return The slot `normalizedTpmMatrix` of the provided [`CAGEexp`] object
#' will be occupied by normalized CAGE signal values per CTSS across all
#' experiments, or with the raw tag counts (in case `method = "none"`).
#' 
#' @references
#' 
#' Balwierz _et al._ (2009) Methods for analyzing deep sequencing expression data:
#' constructing the human and mouse promoterome with deepCAGE data, _Genome Biology_
#' **10**(7):R79. 
#' 
#' @author Vanja Haberle
#' 
#' @seealso [`plotReverseCumulatives`], [`CTSSnormalizedTpmDF`]
#' 
#' @family CAGEr object modifiers
#' @family CAGEr normalised data functions
#' 
#' @examples 
#' ce1 <- normalizeTagCount(exampleCAGEexp, method = "simpleTpm")
#' ce2 <- normalizeTagCount(exampleCAGEexp, method = "powerLaw")
#' 
#' @export

setGeneric( "normalizeTagCount"
          , function( object, method = c("powerLaw", "simpleTpm", "none")
                    , fitInRange = c(10, 1000), alpha = 1.25, T = 10^6)
              standardGeneric("normalizeTagCount"))

#' .normalizeTagCount_switcher
#' 
#' Common code to normalizeTagCount for `CAGEexp` objects
#' Do not reuse elsewhere.
#' 
#' @param method The method.
#' @param object A `CAGEexp` object.
#' 
#' @return A [`DataFrame`].
#' 
#' @noRd

.normalizeTagCount_switcher <- function( method = c("powerLaw", "simpleTpm", "none")
                                       , object, fitInRange, alpha, T) {
  method <- match.arg(method)
  message("\nNormalizing tag count...")
  switch( method
        , powerLaw  = .powerLaw(CTSStagCountDF(object), fitInRange, alpha, T)
        , simpleTpm = .simpleTpm(CTSStagCountDF(object))
        , none      = CTSStagCountDF(object)
        , stop('"method" must be one of ("powerLaw", "simpleTpm", "none")'))
}

# For the CAGEexp class, normalizeTagCount populates the normalized slot or the tagCountMatrix
# experiment.

#' @rdname normalizeTagCount
#' 
setMethod("normalizeTagCount", "CAGEexp", function (object, method, fitInRange, alpha, T) {
	assays(CTSStagCountSE(object), withDimnames=FALSE)$normalizedTpmMatrix <-
	  .normalizeTagCount_switcher(method, object, fitInRange, alpha, T)
	object
})
