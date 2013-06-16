######################################################################################
# Funtions for normalizing CAGE tag count to a referent power-law distribution
# Reference: Balwierz, P. J., Carninci, P., Daub, C. O., Kawai, J., Hayashizaki, 
# Y., Van Belle, W., Beisel, C., et al. (2009). Methods for analyzing deep sequencing 
# expression data: constructing the human and mouse promoterome with deepCAGE data. 
# Genome Biology, 10(7), R79.
#


# function that fits power-law distribution to reverse cumulative of given values - fitting is done using only the range specified in val.range
# ARGUMENTS: values - numerical values whose reverse cumulative distribution will be fitted to power-law (e.g. tag count or signal for regions, peaks, etc.)
#			 val.range - range in which the fitting is done (values outside of this range will not be considered for fitting)
# RETURNS: two coefficients describing fitted power-law distribution (a = slope in the logX-logY plot (where X=signal, Y=number of sites that have >= of given signal), b = intercept in the logX-logY plot)

.fit.power.law.to.reverse.cumulative <- function(values, val.range = c(10, 1000)) {

# using data.table package	
	v <- data.table(num = 1, nr_tags = values)
	v <- v[, sum(num), by = nr_tags]
	setkey(v, nr_tags)
	v$V1 <- rev(cumsum(rev(v$V1)))
	setnames(v, c('nr_tags', 'reverse_cumulative'))
	v <- v[nr_tags >= min(val.range) & nr_tags <= max(val.range)]

#	v <- aggregate(values, by = list(values), FUN = length)
#	v$x <- rev(cumsum(rev(v$x)))
#	colnames(v) <- c('nr_tags', 'reverse_cumulative')
#	v <- subset(v, nr_tags >= min(val.range) & nr_tags <= max(val.range))
	
	lin.m <- lm(log(reverse_cumulative) ~ log(nr_tags), data = v)
	a <- coefficients(lin.m)[2]
	b <- coefficients(lin.m)[1]
	return(c(a, b))
	
}


# function that normalizes values fitted to power-law distribution to a referent power-law distribution
# ARGUMENTS: values - numerical values whose reverse cumulative distribution was fitted to power-law and that will be normalized to referent power-law
#            lin.reg.coef - two coefficients describing power-law distribution fitted to values (as returned by '.fit.power.law.to.reverse.cumulative' function)
#			 alpha - slope of the referent power-law distribution (the actual slope has negative sign and will be -1*alpha)
#			 T - total number of tags (signal) in the referent power-law distribution
# RETURNS: normalized values (vector of the same length as input values); i.e. what would be the value of input values in the referent distribution 

.normalize.to.reference.power.law.distribution <- function(values, lin.reg.coef, alpha = 1.25, T = 10^6) {
	
	a <- lin.reg.coef[1]
	b <- lin.reg.coef[2] 
	lambda <- (T/(zeta(alpha) * exp(b)))^(1/alpha)
	beta <- -1 * a/alpha
	values.norm <- lambda * (values)^beta
	return(values.norm)
	
}
