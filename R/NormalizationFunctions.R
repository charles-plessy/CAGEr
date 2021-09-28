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
  fit.coef <- .fit.power.law.to.reverse.cumulative(values = tag.counts, val.range = fitInRange)
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

#' @name .fit.power.law.to.reverse.cumulative
#' 
#' Function that fits power-law distribution to reverse cumulative of given
#' values - fitting is done using only the range specified in val.range
#' 
#' @return Two coefficients describing fitted power-law distribution (a = slope
#' in the logX-logY plot (where X=signal, Y=number of sites that have >= of given
#' signal), b = intercept in the logX-logY plot)
#' 
#' @importFrom data.table data.table setkeyv setnames
#' @importFrom stats coefficients lm
#' @noRd

.fit.power.law.to.reverse.cumulative <- function(values, val.range = c(10, 1000)) {

# using data.table package	
	v <- data.table(num = 1, nr_tags = values)
	num <- nr_tags <- NULL # Keep R CMD check happy.
	v <- v[, sum(num), by = nr_tags]
	setkeyv(v, "nr_tags")
	v$V1 <- rev(cumsum(rev(v$V1)))
	setnames(v, c('nr_tags', 'reverse_cumulative'))
	
	# check if range values are > 0
	if(any(val.range < 1)){
	  stop(paste("The range of values for fitting the power law",
	             "arg 'fitInRange', expects integers > 0."))
	}
	
	v <- v[nr_tags >= min(val.range) & nr_tags <= max(val.range)]
	
	# check if specified range values have at least 1 entry in v 
	if(nrow(v) < 1){
	  stop(paste("Selected range for fitting the power law does not contain any",
	             "tag count values. Consider changing/increasing 'fitInRange'."))
	}
#	v <- aggregate(values, by = list(values), FUN = length)
#	v$x <- rev(cumsum(rev(v$x)))
#	colnames(v) <- c('nr_tags', 'reverse_cumulative')
#	v <- subset(v, nr_tags >= min(val.range) & nr_tags <= max(val.range))
	
	lin.m <- lm(log(reverse_cumulative) ~ log(nr_tags), data = v)
	
	a <- coefficients(lin.m)[2]
	b <- coefficients(lin.m)[1]
	
	# check if specified range values have >1 entries in v 
	if(is.na(a) && b == 0){
	  stop(paste("Selected range for fitting the power law does not", 
	       "contain enough values. Consider changing/increasing 'fitInRange'."))
	}
	return(c(a, b))
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
