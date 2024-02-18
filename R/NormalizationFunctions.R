######################################################################################
# Functions for normalizing CAGE tag counts

#' .powerLaw
#' 
#' Private funtion for normalizing CAGE tag count to a referent power-law distribution.
#' 
#' @rdname powerLaw
#' 
#' @references
#' Balwierz, P. J., Carninci, P., Daub, C. O., Kawai, J., Hayashizaki,
#' Y., Van Belle, W., Beisel, C., et al. (2009). Methods for analyzing deep sequencing
#' expression data: constructing the human and mouse promoterome with deepCAGE data.
#' Genome Biology, 10(7), R79.
#' 
#' @param tag.counts Numerical values whose reverse cumulative distribution will be fitted
#'        to power-law (e.g. tag count or signal for regions, peaks, etc.)
#' @param fitInRange Range in which the fitting is done (values outside of this range will
#'        not be considered for fitting)
#' @param alpha Slope of the referent power-law distribution (the actual slope has negative
#'        sign and will be -1*alpha)
#' @param T total number of tags (signal) in the referent power-law distribution.
#' 
#' @details S4 Methods are provided for integer vectors, Rle objects, data.frame objects
#' and DataFrame objects, so that the most complex objects can be deconstructed
#' in simpler parts, normalized and reconstructed.
#' 
#' @return Normalized values (vector of the same length as input values); i.e. what would
#' be the value of input values in the referent distribution.  Ouptut objects are numeric,
#' possibly \code{Rle}-encoded or wrapped in \code{data.frames} or \code{DataFrames}
#' according to the input.

setGeneric(".powerLaw", function(tag.counts, fitInRange = c(10, 1000), alpha = 1.25, T = 10^6) {
  standardGeneric(".powerLaw")})

setMethod(".powerLaw", "numeric", function (tag.counts, fitInRange, alpha, T) {
  fit.coef <- fitPowerLaw(tag.counts, fitInRange)
	.normalize.to.reference.power.law.distribution(values = tag.counts, lin.reg.coef = fit.coef, alpha = alpha, T = T)
})

setMethod(".powerLaw", "Rle", function (tag.counts, fitInRange, alpha, T) {
  Rle(.powerLaw(as.integer(tag.counts), fitInRange, alpha, T))
})

setMethod(".powerLaw", "data.frame", function (tag.counts, fitInRange, alpha, T) {
  data.frame(sapply(tag.counts, function(X) .powerLaw(X, fitInRange, alpha, T)))
})

setMethod(".powerLaw", "DataFrame", function (tag.counts, fitInRange, alpha, T) {
  DataFrame(sapply(tag.counts, function(X) .powerLaw(X, fitInRange, alpha, T)))
})

#' Fit a power law to the CAGE data
#' 
#' Function that fits power-law distribution to reverse cumulative of given
#' values.  Fitting is done using only a specified range.
#' 
#' @param x The object containing the data.
#' @param fitInRange The range in which to fit the power law.
#' 
#' @return Returns the slope (`plSlope`) and the intercept (`plInt`) of the
#' fitted power-law distribution.  If the input is a single sample, function
#' returns a numeric vector.  If the input is multiple samples, the function
#' returns a `DataFrame`.  If the input is a whole _CAGEexp_ object, the
#' function returns the object with information added to its `colData`.
#' 
#' @importFrom stats coefficients lm
#' 
#' @author Vanja Haberle
#' @author Charles Plessy
#' 
#' @export

setGeneric("fitPowerLaw", function(x, fitInRange = c(10, 1000)) {
  standardGeneric("fitPowerLaw")})

#' @rdname fitPowerLaw
#' @examples
#' exampleCAGEexp |> fitPowerLaw() |> colData()

setMethod("fitPowerLaw", "CAGEexp", function (x, fitInRange) {
  DF <- fitPowerLaw(CTSStagCountDF(x), fitInRange)
  colData(x) <- cbind(colData(x), DF)
  x
})

#' @rdname fitPowerLaw
#' @examples
#' exampleCAGEexp |> CTSStagCountDF() |> fitPowerLaw()

setMethod("fitPowerLaw", "DataFrame", function (x, fitInRange) {
  m <- sapply(x, fitPowerLaw, fitInRange)
  DF <- as(t(m), "DataFrame")
  rownames(DF) <- names(x)
  DF
})

#' @rdname fitPowerLaw
#' @examples
#' exampleCAGEexp |> CTSStagCountGR('all') |> fitPowerLaw()

setMethod("fitPowerLaw", "GRangesList", function (x, fitInRange) {
  m <- sapply(x, fitPowerLaw, fitInRange)
  DF <- as(t(m), "DataFrame")
  rownames(DF) <- names(x)
  DF
})

#' @rdname fitPowerLaw
#' @examples
#' exampleCAGEexp |> CTSStagCountGR(1) |> fitPowerLaw()

setMethod("fitPowerLaw", "GRanges", function (x, fitInRange) {
  fitPowerLaw(score(x), fitInRange)
})

#' @rdname fitPowerLaw
#' @examples
#' exampleCAGEexp |> CTSStagCountGR(1) |> score() |> fitPowerLaw()

setMethod("fitPowerLaw", "Rle", function (x, fitInRange) {
  .fit.power.law.to.reverse.cumulative(x, fitInRange)
})

#' @rdname fitPowerLaw
#' @examples
#' exampleCAGEexp |> CTSStagCountGR(1) |> score() |> decode() |> fitPowerLaw()

setMethod("fitPowerLaw", "numeric", function(x, fitInRange) {
  .fit.power.law.to.reverse.cumulative(x, fitInRange)
})

.fit.power.law.to.reverse.cumulative <- function(values, val.range = c(10, 1000)) {
  # Using Rle
  v <- sort(Rle(values), decr = TRUE)
  nr_tags <- rev(runValue(v))
  reverse_cumulative <- v  |> runLength() |> cumsum() |> rev()
  
  inRange <- nr_tags >= min(val.range) & nr_tags <= max(val.range)
  
  lin.m <- lm(log(reverse_cumulative[inRange]) ~ log(nr_tags[inRange]))
  
  a <- coefficients(lin.m)[2]
  names(a) <- 'plSlope'
  b <- coefficients(lin.m)[1]
  names(b) <- 'plInt'
  
  # check if specified range values have >1 entries in v 
  if(is.na(a) && b == 0){
    stop(paste("Selected range for fitting the power law does not", 
               "contain enough values. Consider changing/increasing 'fitInRange'."))
  }
  c(a, b)
}

#' .normalize.to.reference.power.law.distribution
#' 
#' Function that normalizes values fitted to power-law distribution to a referent
#' power-law distribution
#' 
#' @param values - numerical values whose reverse cumulative distribution was fitted to
#'        power-law and that will be normalized to referent power-law
#' @param lin.reg.coef - two coefficients describing power-law distribution fitted to
#'        values (as returned by '.fit.power.law.to.reverse.cumulative' function)
#' 
#' @importFrom VGAM zeta
#' @noRd

.normalize.to.reference.power.law.distribution <- function(values, lin.reg.coef, alpha = 1.25, T = 10^6) {
	
	a <- lin.reg.coef[1]
	b <- lin.reg.coef[2] 
	lambda <- (T/(VGAM::zeta(alpha) * exp(b)))^(1/alpha)
	beta <- -1 * a/alpha
	values.norm <- lambda * (values)^beta
	return(values.norm)
	
}

#' @name .simpleTpm
#' @noRd
#' 
#' @param tag.counts An object containing tag counts.
#' 
#' Methods are provided for integer vectors, Rle objects, data.frame objects
#' and DataFrame objects, so that the most complex objects can be deconstructed
#' in simpler parts, normalized and reconstructed.
#' 
#' @return Integers become numeric, Rle, data.frame and DataFrame are conserved.

setGeneric(".simpleTpm", function(tag.counts) {
  standardGeneric(".simpleTpm")})

setMethod(".simpleTpm", "numeric", function (tag.counts) {
  tag.counts / sum(tag.counts) * 10^6
})

setMethod(".simpleTpm", "Rle", function (tag.counts) {
  Rle(.simpleTpm(as.integer(tag.counts)))
})

setMethod(".simpleTpm", "data.frame", function (tag.counts) {
  data.frame(sapply(tag.counts, .simpleTpm))
})

setMethod(".simpleTpm", "DataFrame", function (tag.counts) {
  DataFrame(sapply(tag.counts, .simpleTpm))
})
