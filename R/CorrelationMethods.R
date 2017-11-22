#' @include CAGEr.R

#' @name plotCorrelation
#' 
#' @title Pairwise scatter plots and correlations of CAGE signal
#' 
#' @description Calculates the pairwise correlation between samples and creates
#' a plot matrix showing the correlation coeficients in the upper triangle, the
#' sample names in the diagonal, and the catter plots in the lower triangle.
#' 
#' @param object A \code{\link{CAGEr}} object.
#' 
#' @param what The clustering level to be used for plotting and calculating
#'   correlations.  Can be either \code{"CTSS"} to use individual TSSs or
#'   \code{"consensusClusters"} to use consensus clusters, \emph{i.e.} entire
#'   promoters.
#' 
#' @param values Use either \code{"raw"} (default) or \code{"normalized"} CAGE
#'   signal.
#' 
#' @param samples Character vector indicating which samples to use.  Can be
#'   either \code{"all"} to select all samples in a \code{CAGEr} object, or a
#'   subset of valid sample labels as returned by the
#'   \code{\link{sampleLabels}} function.
#' 
#' @param method A character string indicating which correlation coefficient
#'   should be computed.  Passed to \code{cor} function.  Can be one of
#'   \code{"pearson"}, \code{"spearman"}, or \code{"kendall"}.
#' 
#' @param tagCountThreshold Only TSSs with tag count \code{>= tagCountThreshold}
#'   in either one (\code{applyThresholdBoth = FALSE}) or both samples
#'   (\code{applyThresholdBoth = TRUE}) are plotted and used to calculate
#'   correlation.
#' 
#' @param applyThresholdBoth See \code{tagCountThreshold} above.
#' 
#' @param plotSize Size of the individual comparison plot in pixels - the
#' total size of the resulting png will be \code{length(samples) * plotSize}
#' in both dimensions.  Ignored in \code{plotCorrelation2}.
#' 
#' @details Internally, the data is converted to a \code{data.frame} format,
#' so there may be performance issues if there are many samples.  On the other
#' hand, this kind of plot does not make much sense for large numbers of samples.
#' 
#' In the scatter plots, a pseodo-count equal to half the lowest score is added
#' to the null values so that they can appear despite logarithmic scale.
#' 
#' @return Displays the plot and returns a \code{matrix} of pairwise
#' correlations between selected samples.  The scatterplots of
#' \code{plotCorrelation} are colored according to the density of points, and
#' in \code{plotCorrelation2} they are just black and white, which is much
#' faster to plot.  Note that while the scatterplots are on a logarithmic scale
#' with pseudocount added to the zero values, the correlation coefficients are
#' calculated on untransformed (but thresholded) data.
#' 
#' @author Vanja Haberle
#' @author Charles Plessy
#' 
#' @family CAGEr plot functions
#' 
#' @importFrom graphics axis
#' @importFrom S4Vectors decode
#' 
#' @examples
#' plotCorrelation(object = exampleCAGEset)
#' 
#' plotCorrelation2(exampleCAGEexp, what = "consensusClusters", value = "normalized")
#' 
#' @export

setGeneric( "plotCorrelation"
          , function( object, what = c("CTSS", "consensusClusters")
                    , values = c("raw", "normalized")
                    , samples = "all", method = "pearson"
                    , tagCountThreshold = 1, applyThresholdBoth = FALSE, plotSize=800)
              standardGeneric("plotCorrelation"))

#' @rdname plotCorrelation

setMethod( "plotCorrelation", "CAGEr"
         , function( object, what, values, samples, method
                   , tagCountThreshold, applyThresholdBoth, plotSize) {
  what   <- match.arg(what)
  values <- match.arg(values)
	if (what == "CTSS" & values == "raw")
		tag.count <- CTSStagCountDF(object)
  if (what == "CTSS" & values == "normalized")
		tag.count <- CTSSnormalizedTpmDF(object)
	if (what == "consensusClusters" & values == "raw")
		stop("Raw consensus clusters not supported yet.")
	if (what == "consensusClusters" & values == "normalized")
		tag.count <- consensusClustersTpm(object)

	if(all(samples %in% sampleLabels(object))){
		tag.count <- tag.count[,samples]
		nr.samples <- length(samples)
	}else if(samples == "all"){
		samples <- sampleLabels(object)
		nr.samples <- length(samples)
	}else{
		stop("'samples' parameter must be either \"all\" or a character vector of valid sample labels!")
	}
	
	corr.m <- matrix(rep(1, (nr.samples)^2), nrow = nr.samples)
	colnames(corr.m) <- samples
	rownames(corr.m) <- samples
	
	old.par <- par(mfrow = c(nr.samples, nr.samples), mai = c(0.05,0.05,0.05,0.05), omi = c(0.05*plotSize*(log2(nr.samples)/log2(3) + (2-1/log2(3)))/360,0.1*plotSize*(log2(nr.samples)/log2(3) + (2-1/log2(3)))/360,0,0))
	on.exit(par(old.par))
	
	for(i in c(1:nr.samples)){
		for(j in c(1:nr.samples)){
			
			if(i == j){
				
				plot(1, 1, type = "n", bty = "n", xlim = c(0,1), ylim = c(0,1), axes = F)
				text(0.5, 0.5, samples[i], cex = 0.98/max(sapply(samples, strwidth)))
				box(lwd = 3)
		
			}else {
			
				x <- tag.count[,samples[j]]
				y <- tag.count[,samples[i]]
				if(applyThresholdBoth){
					idx <- (x >= tagCountThreshold) & (y >= tagCountThreshold)
				}else{
					idx <- (x >= tagCountThreshold) | (y >= tagCountThreshold)
				}
				x <- x[idx]
				y <- y[idx]
			
				if (j > i) {
					pairwise.cor <- cor(x = decode(x), y = decode(y), method = method) # decode() in case of Rle.
					plot(1, 1, type = "n", bty = "n", xlim = c(0,1), ylim = c(0,1), axes = F)
					txt <- sprintf("%.2f", pairwise.cor)
					txt.abs <- sprintf("%.2f", abs(pairwise.cor))
					text(0.5, 0.5, txt, cex = 1.5 + 0.5/strwidth(txt.abs) * abs(pairwise.cor))
					box(lwd = 3)
					corr.m[i,j] <- pairwise.cor
					corr.m[j,i] <- pairwise.cor
			
				}else{
					.mySmoothScatter(x = log10(x+1), y = log10(y+1), xlim = c(0, 3), ylim = c(0,3), nrpoints = 0, nbin = c(plotSize, plotSize), bandwidth = c(3/plotSize * 5, 3/plotSize * 5), transformation = function(x) x^(1/6), axes = F)
					if(i == nr.samples & j < nr.samples){
						if((nr.samples <= 3) | ((nr.samples > 3) & (j%%2 == 1))){
							axis(side = 1, at = seq(0,3), labels = 10^seq(0,3), cex.axis = 0.1*plotSize*(log2(nr.samples)/log2(3) + (2-1/log2(3)))/360/(4*strwidth("1")))
						}else{
							axis(side = 1, at = seq(0,3), labels = rep("", 4))
						}
					}
					if(j == 1 & i > 1){
						if((nr.samples <= 3) | ((nr.samples > 3) & ((i - nr.samples)%%2 == 0))){
							axis(side = 2, at = seq(0,3), labels = 10^seq(0,3), las = 2, cex.axis = 0.1*plotSize*(log2(nr.samples)/log2(3) + (2-1/log2(3)))/360/(4*strwidth("1")))
						}else{
							axis(side = 2, at = seq(0,3), labels = rep("", 4))
						}
					}
					box(lwd = 3)
				}
#				fit.lm <- lm(y ~ x)
#				.mySmoothScatter(x = x, y = y, xlim = c(0, 200), ylim = c(0,200), nrpoints = 0, nbin = c(800, 800), bandwidth = c(3, 3), transformation = function(x) x^(1/9), axes = F)
#				lines(x = c(0,10,100,1000), y = coefficients(fit.lm)[2]*c(0,10,100,1000) + coefficients(fit.lm)[1], col = "red3", lwd = 3)				
				
			}
		}
	}
	
	return(corr.m)
})

#' @rdname plotCorrelation
#' @export

setGeneric( "plotCorrelation2"
          , function( object, what = c("CTSS", "consensusClusters")
                    , values = c("raw", "normalized")
                    , samples = "all", method = "pearson"
                    , tagCountThreshold = 1, applyThresholdBoth = FALSE)
              standardGeneric("plotCorrelation2"))

#' @importFrom graphics pairs
#' @rdname plotCorrelation

setMethod( "plotCorrelation2", "CAGEr"
         , function( object, what, values, samples, method
                   , tagCountThreshold, applyThresholdBoth) {
  what   <- match.arg(what)
  values <- match.arg(values)
  
  # Select input data
  if (what == "CTSS" & values == "raw")
    tag.count <- CTSStagCountDF(object)
  if (what == "CTSS" & values == "normalized")
    tag.count <- CTSSnormalizedTpmDF(object)
  if (what == "consensusClusters" & values == "raw")
    tag.count <- assay(consensusClustersSE(object), "counts")
  if (what == "consensusClusters" & values == "normalized")
    tag.count <- assay(consensusClustersSE(object), "normalized")
  
  # Coerce DataFrame to matrix
  if (class(tag.count) == "DataFrame")
    tag.count <- as.matrix(as.data.frame(tag.count))
  
  # Select samples
  if (all(samples %in% sampleLabels(object))) {
    tag.count <- tag.count[,samples]
  } else if(samples == "all"){
    samples <- sampleLabels(object)
  } else stop("'samples' parameter must be either \"all\" or a character vector of valid sample labels!")
  nr.samples <- length(samples)
  
  # Function to apply threshold pairwise
  applyThreshold <- function(df) {
    if (applyThresholdBoth) {
      idx <- (df[,1] >= tagCountThreshold) & (df[,2] >= tagCountThreshold)
    } else {
      idx <- (df[,1] >= tagCountThreshold) | (df[,2] >= tagCountThreshold)
    }
    df[idx,]
  }
  
  # Pre-calculate a vector of correlation coefficients
  corTreshold <- function(x,y) {
    df <- data.frame(x, y)
    df <- applyThreshold(df)
    cor(x = df$x, y = df$y, method = method)
  }
  
  corr.v <- numeric()
  for (i in 1:(nr.samples-1)) {
    for (j in (min(i+1, nr.samples)):nr.samples) {
      corr.v <- append(corr.v, corTreshold(tag.count[,i], tag.count[,j]))
    }
  }
  
  # Add pseudocount to null values so that the plot axes are correctly set.
  nulls <- tag.count == 0
  pseudocount      <- min(tag.count[! nulls]) / 2
  tag.count[nulls] <- tag.count[nulls] + pseudocount
  
  # This closure retreives correlation coefficients one after the other.
  mkPanelCor <- function() {
    i <- 1
    function(x, y, digits=2, prefix="", cex.cor, ...) {
      r <- corr.v[i]
      i <<- i + 1
      txt <- format(c(r, 0.123456789), digits=digits)[1]
      txt <- paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      gmean <- function(x) {
        x <- x[x > 0]
        exp(mean(log(range(x))))
      }
      text(gmean(x), gmean(y), txt, cex = cex.cor * sqrt(r))
    }
  }
  panel.cor <- mkPanelCor()

  pointsUnique <- function(x,y,...) {
    df <- data.frame(x, y)
    df <- unique(df)
    df <- applyThreshold(df)
    df <- df[rowSums(df) > pseudocount * 2,] # Remove the (0,0) point
    points(df, ...)
  }

  pairs( tag.count
       , lower.panel = pointsUnique
       , upper.panel = panel.cor
       , pch = "."
       , cex = 4
       , log = "xy"
       , las = 1
       , xaxp = c(1,10,1)
       , yaxp = c(1,10,1))
  
  # Return a correlation matrix
  corr.m <- matrix(0, nr.samples, nr.samples)
  colnames(corr.m) <- samples
  rownames(corr.m) <- samples
  corr.m[upper.tri(corr.m)] <- corr.v
  corr.m[lower.tri(corr.m)] <- t(corr.m)[lower.tri(corr.m)]
  corr.m
})


# my version of smooth scatter that allows passing range.x argument to grDevices:::.smoothScatterCalcDensity function to calculate 2D kernel smoothed density

#' @importFrom grDevices blues9 colorRamp colorRampPalette xy.coords
#' @importFrom graphics points
#' @importFrom KernSmooth bkde2D

# Local copy of grDevices:::.smoothScatterCalcDensity,
# to avoid problems in case the original function is changed
# (since the original is private, we can not assume that changes maintain
#  compatibility with existing code.)
grDevices.smoothScatterCalcDensity <- function (x, nbin, bandwidth, range.x) 
{
    if (length(nbin) == 1) 
        nbin <- c(nbin, nbin)
    if (!is.numeric(nbin) || length(nbin) != 2) 
        stop("'nbin' must be numeric of length 1 or 2")
    if (missing(bandwidth)) {
        bandwidth <- diff(apply(x, 2, stats::quantile, probs = c(0.05, 
            0.95), na.rm = TRUE, names = FALSE))/25
        bandwidth[bandwidth == 0] <- 1
    }
    else {
        if (!is.numeric(bandwidth)) 
            stop("'bandwidth' must be numeric")
        if (any(bandwidth <= 0)) 
            stop("'bandwidth' must be positive")
    }
    rv <- KernSmooth::bkde2D(x, bandwidth = bandwidth, gridsize = nbin, 
        range.x = range.x)
    rv$bandwidth <- bandwidth
    rv
}

.mySmoothScatter <- function (x, y = NULL, nbin = 128, bandwidth, colramp = colorRampPalette(c("white",
blues9)), nrpoints = 100, pch = ".", cex = 1, col = "black",
transformation = function(x) x^0.25, postPlotHook = box,
xlab = NULL, ylab = NULL, xlim, ylim, xaxs = par("xaxs"),
yaxs = par("yaxs"), ...)
{
    if (!is.numeric(nrpoints) | (nrpoints < 0) | (length(nrpoints) !=
    1))
    stop("'nrpoints' should be numeric scalar with value >= 0.")
    xlabel <- if (!missing(x))
    deparse(substitute(x))
    ylabel <- if (!missing(y))
    deparse(substitute(y))
    xy <- xy.coords(x, y, xlabel, ylabel)
    xlab <- if (is.null(xlab))
    xy$xlab
    else xlab
    ylab <- if (is.null(ylab))
    xy$ylab
    else ylab
    x <- cbind(xy$x, xy$y)[is.finite(xy$x) & is.finite(xy$y),
    , drop = FALSE]
    if (!missing(xlim)) {
        stopifnot(is.numeric(xlim), length(xlim) == 2, is.finite(xlim))
        x <- x[min(xlim) <= x[, 1] & x[, 1] <= max(xlim), ]
    }
    else {
        xlim <- range(x[, 1])
    }
    if (!missing(ylim)) {
        stopifnot(is.numeric(ylim), length(ylim) == 2, is.finite(ylim))
        x <- x[min(ylim) <= x[, 2] & x[, 2] <= max(ylim), ]
    }
    else {
        ylim <- range(x[, 2])
    }
    map <- grDevices.smoothScatterCalcDensity(x, nbin, bandwidth, range.x = list(xlim = c(xlim[1] - 1.5*bandwidth[1], xlim[2] + 1.5*bandwidth[1]), ylim = c(ylim[1] - 1.5*bandwidth[2], ylim[2] + 1.5*bandwidth[2])))
    xm <- map$x1
    ym <- map$x2
    dens <- map$fhat
    dens[] <- transformation(dens)
    image(xm, ym, z = dens, col = colramp(256), xlab = xlab,
    ylab = ylab, xlim = xlim, ylim = ylim, xaxs = xaxs, yaxs = yaxs,
    ...)
    if (!is.null(postPlotHook))
    postPlotHook()
    if (nrpoints > 0) {
        nrpoints <- min(nrow(x), ceiling(nrpoints))
        stopifnot((nx <- length(xm)) == nrow(dens), (ny <- length(ym)) ==
        ncol(dens))
        ixm <- 1L + as.integer((nx - 1) * (x[, 1] - xm[1])/(xm[nx] -
        xm[1]))
        iym <- 1L + as.integer((ny - 1) * (x[, 2] - ym[1])/(ym[ny] -
        ym[1]))
        sel <- order(dens[cbind(ixm, iym)])[seq_len(nrpoints)]
        points(x[sel, ], pch = pch, cex = cex, col = col)
    }
}
