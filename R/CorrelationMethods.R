#' @include CAGEr.R

#' @name plotCorrelation
#' 
#' @title Pairwise scatter plots and correlations of CAGE signal
#' 
#' @description Calculates the pairwise correlation between samples and creates
#' a plot matrix showing the correlation coeficients in the upper triangle, the
#' sample names in the diagonal, and the catter plots in the lower triangle.
#' 
#' @param object A \code{\link{CAGEr}} object or (only for
#'   \code{plotCorrelation2}) a \code{\link{SummarizedExperiment}} or an
#'   expression table as a \code{\link{DataFrame}}, \code{\link{data.frame}} or
#'   \code{\link{matrix}} object.
#' 
#' @param what The clustering level to be used for plotting and calculating
#'   correlations.  Can be either \code{"CTSS"} to use individual TSSs or
#'   \code{"consensusClusters"} to use consensus clusters, \emph{i.e.} entire
#'   promoters.  Ignored for anything else than \code{CAGEr} objects.
#' 
#' @param values Use either \code{"raw"} (default) or \code{"normalized"} CAGE
#'   signal.  Ignored for plain expression tables.
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
#' @details In the scatter plots, a pseudo-count equal to half the lowest score
#' is added to the null values so that they can appear despite logarithmic scale.
#' 
#' \code{SummarizedExperiment} objects are expected to contain raw tag counts
#' in a \dQuote{counts} assay and the normalized expression scores in a
#' \dQuote{normalized} assay.
#' 
#' Avoid using large \code{matrix} objects as they are coerced to
#' \code{DataFrame} class without special care for efficiency.
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
#' 
#' @param digits The number of significant digits for the data to be kept in log
#'   scale.  Ignored in \code{plotCorrelation}.  In \code{plotCorrelation2}, the
#'   number of points plotted is considerably reduced by rounding the point
#'   coordinates to a small number of significant digits before removing
#'   duplicates.  Chose a value that makes the plot visually indistinguishable
#'   with non-deduplicated data, by making tests on a subset of the data.
#' 
#' @details \code{plotCorrelation2} speeds up the plotting by a) deduplicating
#' that data: no point is plot twice at the same coordinates, b) rounding the
#' data so that indistinguishable positions are plotted only once, c) using a
#' black square glyph for the points, d) caching some calculations that are
#' made repeatedly (to determine where to plot the correlation coefficients),
#' and e) preventing coercion of \code{DataFrames} to \code{data.frames}.
#' 
#' @importFrom memoise memoise
#' @export

setGeneric( "plotCorrelation2"
          , function( object, what = c("CTSS", "consensusClusters")
                    , values = c("raw", "normalized")
                    , samples = "all", method = "pearson"
                    , tagCountThreshold = 1, applyThresholdBoth = FALSE
                    , digits = 3)
              standardGeneric("plotCorrelation2"))

#' @importFrom graphics pairs
#' @rdname plotCorrelation

setMethod( "plotCorrelation2", "CAGEexp"
         , function( object, what, values, samples, method
                   , tagCountThreshold, applyThresholdBoth
                   , digits) {
  what <- match.arg(what)
  se <- switch( what
              , CTSS              = CTSStagCountSE(object)
              , consensusClusters = consensusClustersSE(object)
              , genes             = GeneExpSE(object)
              , stop("Unsupported value for ", dQuote("what"), ": ", what))
  plotCorrelation2( se
                  , what               = what
                  , values             = values
                  , samples            = samples
                  , method             = method
                  , tagCountThreshold  = tagCountThreshold
                  , applyThresholdBoth = applyThresholdBoth
                  , digits             = digits)
})

#' @rdname plotCorrelation

setMethod( "plotCorrelation2", "SummarizedExperiment"
         , function( object, what, values, samples, method
                   , tagCountThreshold, applyThresholdBoth
                   , digits) {
  values <- match.arg(values)
  
  if (values == "raw") {
    if ("counts" %in% assayNames(object)) {
      values <- "counts"
    } else {
      stop( "Could not find a ", dQuote("counts"), " assay for the "
          , dQuote(what), " clustering level")
    }
  } else if (values == "normalized") {
    if ("normalized" %in% assayNames(object)) {
      values <- "normalized"
    } else if ("normalizedTpmMatrix" %in% assayNames(object)) {
      values <- "normalizedTpmMatrix"
    } else {
      stop( "Could not find a ", dQuote("normalized"), " assay for the "
          , dQuote(what), " clustering level")
    }
  }
  
  plotCorrelation2( assay(object, values)
                  , what               = what
                  , values             = values
                  , samples            = samples
                  , method             = method
                  , tagCountThreshold  = tagCountThreshold
                  , applyThresholdBoth = applyThresholdBoth
                  , digits             = digits)
})

#' @rdname plotCorrelation

setMethod( "plotCorrelation2", "DataFrame"
         , function( object, what, values, samples, method
                   , tagCountThreshold, applyThresholdBoth
                   , digits) {
  .plotCorrelation2( object
                   , samples            = samples
                   , method             = method
                   , tagCountThreshold  = tagCountThreshold
                   , applyThresholdBoth = applyThresholdBoth
                   , digits             = digits)
})

#' @rdname plotCorrelation

setMethod( "plotCorrelation2", "data.frame"
         , function( object, what, values, samples, method
                   , tagCountThreshold, applyThresholdBoth
                   , digits) {
  .plotCorrelation2( object
                   , samples            = samples
                   , method             = method
                   , tagCountThreshold  = tagCountThreshold
                   , applyThresholdBoth = applyThresholdBoth
                   , digits             = digits)
})

#' @rdname plotCorrelation

setMethod( "plotCorrelation2", "matrix"
         , function( object, what, values, samples, method
                   , tagCountThreshold, applyThresholdBoth
                   , digits) {
  .plotCorrelation2( as.data.frame(object)
                   , samples            = samples
                   , method             = method
                   , tagCountThreshold  = tagCountThreshold
                   , applyThresholdBoth = applyThresholdBoth
                   , digits             = digits)
})


# Helper function to apply threshold pairwise
.applyThreshold <- function(df, tagCountThreshold, applyThresholdBoth) {
  if (applyThresholdBoth) {
    idx <- (df[[1]] >= tagCountThreshold) & (df[[2]] >= tagCountThreshold)
  } else {
    idx <- (df[[1]] >= tagCountThreshold) | (df[[2]] >= tagCountThreshold)
  }
  df[idx,]
}

# Helper function to Pre-calculate a vector of correlation coefficients
corVector <- function(expr.table, method, tagCountThreshold, applyThresholdBoth) {
  corTreshold <- function(x, y, method) {
    df <- data.frame(x, y)
    df <- .applyThreshold(df, tagCountThreshold, applyThresholdBoth)
    cor(x = df$x, y = df$y, method = method)
  }
  nr.samples <- ncol(expr.table)
  corr.v <- numeric()
  for (i in 1:(nr.samples-1)) {
    for (j in (min(i+1, nr.samples)):nr.samples) {
      corr.v <- append(corr.v, corTreshold(expr.table[[i]], expr.table[[j]], method))
    }
  }
  corr.v
}

#' @importFrom grDevices dev.flush dev.hold
#' @importFrom graphics Axis mtext

# Re-implement the pairs function to prevent coercion to data.frame
pairs.DataFrame <- function (x, labels, panel = points, ..., horInd = 1:nc, verInd = 1:nc, 
    lower.panel = panel, upper.panel = panel, diag.panel = NULL, 
    text.panel = textPanel, label.pos = 0.5 + has.diag/3, line.main = 3, 
    cex.labels = NULL, font.labels = 1, row1attop = TRUE, gap = 1, 
    log = "") 
{
    if (doText <- missing(text.panel) || is.function(text.panel)) 
        textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
            y, txt, cex = cex, font = font)
    localAxis <- function(side, x, y, xpd, bg, col = NULL, main, 
        oma, ...) {
        xpd <- NA
        if (side%%2L == 1L && xl[j]) 
            xpd <- FALSE
        if (side%%2L == 0L && yl[i]) 
            xpd <- FALSE
        if (side%%2L == 1L) 
            Axis(x, side = side, xpd = xpd, ...)
        else Axis(y, side = side, xpd = xpd, ...)
    }
    localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
    localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
    localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
    localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
    dots <- list(...)
    nmdots <- names(dots)
    # if (!is.matrix(x)) {
    #     x <- as.data.frame(x)
    #     for (i in seq_along(names(x))) {
    #         if (is.factor(x[[i]]) || is.logical(x[[i]])) 
    #             x[[i]] <- as.numeric(x[[i]])
    #         if (!is.numeric(unclass(x[[i]]))) 
    #             stop("non-numeric argument to 'pairs'")
    #     }
    # }
    # else if (!is.numeric(x)) 
    #     stop("non-numeric argument to 'pairs'")
    panel <- match.fun(panel)
    if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
        lower.panel <- match.fun(lower.panel)
    if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
        upper.panel <- match.fun(upper.panel)
    if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) 
        diag.panel <- match.fun(diag.panel)
    if (row1attop) {
        tmp <- lower.panel
        lower.panel <- upper.panel
        upper.panel <- tmp
        tmp <- has.lower
        has.lower <- has.upper
        has.upper <- tmp
    }
    nc <- ncol(x)
    if (nc < 2L) 
        stop("only one column in the argument to 'pairs'")
    if (!all(horInd >= 1L & horInd <= nc)) 
        stop("invalid argument 'horInd'")
    if (!all(verInd >= 1L & verInd <= nc)) 
        stop("invalid argument 'verInd'")
    if (doText) {
        if (missing(labels)) {
            labels <- colnames(x)
            if (is.null(labels)) 
                labels <- paste("var", 1L:nc)
        }
        else if (is.null(labels)) 
            doText <- FALSE
    }
    oma <- if ("oma" %in% nmdots) 
        dots$oma
    main <- if ("main" %in% nmdots) 
        dots$main
    if (is.null(oma)) 
        oma <- c(4, 4, if (!is.null(main)) 6 else 4, 4)
    opar <- par(mfrow = c(length(horInd), length(verInd)), mar = rep.int(gap/2, 
        4), oma = oma)
    on.exit(par(opar))
    dev.hold()
    on.exit(dev.flush(), add = TRUE)
    xl <- yl <- logical(nc)
    if (is.numeric(log)) 
        xl[log] <- yl[log] <- TRUE
    else {
        xl[] <- grepl("x", log)
        yl[] <- grepl("y", log)
    }
    for (i in if (row1attop) 
        verInd
    else rev(verInd)) for (j in horInd) {
        l <- paste0(ifelse(xl[j], "x", ""), ifelse(yl[i], "y", 
            ""))
        localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
            type = "n", ..., log = l)
        if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
            box()
            if (i == 1 && (!(j%%2L) || !has.upper || !has.lower)) 
                localAxis(1L + 2L * row1attop, x[, j], x[, i], 
                  ...)
            if (i == nc && (j%%2L || !has.upper || !has.lower)) 
                localAxis(3L - 2L * row1attop, x[, j], x[, i], 
                  ...)
            if (j == 1 && (!(i%%2L) || !has.upper || !has.lower)) 
                localAxis(2L, x[, j], x[, i], ...)
            if (j == nc && (i%%2L || !has.upper || !has.lower)) 
                localAxis(4L, x[, j], x[, i], ...)
            mfg <- par("mfg")
            if (i == j) {
                if (has.diag) 
                  localDiagPanel(x[, i], ...)
                if (doText) {
                  par(usr = c(0, 1, 0, 1))
                  if (is.null(cex.labels)) {
                    l.wid <- strwidth(labels, "user")
                    cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
                  }
                  xlp <- if (xl[i]) 
                    10^0.5
                  else 0.5
                  ylp <- if (yl[j]) 
                    10^label.pos
                  else label.pos
                  text.panel(xlp, ylp, labels[i], cex = cex.labels, 
                    font = font.labels)
                }
            }
            else if (i < j) 
                localLowerPanel(x[, j], x[, i], ...)
            else localUpperPanel(x[, j], x[, i], ...)
            if (any(par("mfg") != mfg)) 
                stop("the 'panel' function made a new plot")
        }
        else par(new = FALSE)
    }
    if (!is.null(main)) {
        font.main <- if ("font.main" %in% nmdots) 
            dots$font.main
        else par("font.main")
        cex.main <- if ("cex.main" %in% nmdots) 
            dots$cex.main
        else par("cex.main")
        mtext(main, 3, line.main, outer = TRUE, at = 0.5, cex = cex.main, 
            font = font.main)
    }
    invisible(NULL)
}

# The function that runs the actual work of calculating correlations and
# plotting expression values.
.plotCorrelation2 <- function( expr.table, samples, method
                             , tagCountThreshold, applyThresholdBoth
                             , digits) {
  # Select samples
  if (all(samples %in% colnames(expr.table))) {
    expr.table <- expr.table[,samples]
  } else if(samples == "all"){
    samples <- colnames(expr.table)
  } else stop("'samples' parameter must be either \"all\" or a character vector of valid sample labels!")
  nr.samples <- length(samples)
  
  # Pre-calculate a vector of correlation coefficients
  corr.v <- corVector(expr.table, method, tagCountThreshold, applyThresholdBoth)
  
  # Add pseudocount to null values so that the plot axes are correctly set.
  pseudocount <- min(sapply(expr.table, function(x) min(x[x>0]))) / 2
  expr.table  <- DataFrame(lapply( expr.table
                                 , function(x) {x[x==0] <- pseudocount ; x}))
  
  # This closure retreives correlation coefficients one after the other.
  mkPanelCor <- function() {
    i <- 1
    function(x, y, digits=2, prefix="", cex.cor, ...) {
      r <- corr.v[i]
      i <<- i + 1
      txt <- format(c(r, 0.123456789), digits=digits)[1]
      txt <- paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      gmean <- memoise(function(x) {
        exp(mean(log(c(pseudocount, max(x)))))
      })
      text(gmean(x), gmean(y), txt, cex = cex.cor * sqrt(r))
    }
  }
  panel.cor <- mkPanelCor()
  
  # uniqueSignif returns a plain data.frame, because compression is already
  # maximal in this context.  See benchmarking of alternatives algorithms
  # in the file "benchmarks/unique-signif.md" in the CAGEr's Git repository.
  uniqueSignif <- function(x, y, digits = 0, log = c("", "xy")) {
    log <- match.arg(log)
    if (log == "xy") {x <- log1p(x) ; y <- log1p(y)}
    u  <- unique(Rle(complex( real      = decode(signif(x, digits = digits))
                            , imaginary = decode(signif(y, digits = digits)))))
    df <- data.frame(x = Re(u), y = Im(u))
    if (log == "xy") df <- expm1(df)
    df
  }

  # Thresholds are lowered of a minute amont because the rounding in the log
  # scale in uniqueSignif adds minute errors to values that were already round.
  pointsUnique <- function(x,y,...) {
    df <- uniqueSignif(x, y, digits = digits, log = "xy")
    df <- .applyThreshold(df, tagCountThreshold * 0.999, applyThresholdBoth)
    df <- df[rowSums(df) > pseudocount * 1.999,] # Remove the (0,0) point.
    points(df, ...)
  }

  pairs( expr.table
       , lower.panel = pointsUnique
       , upper.panel = panel.cor
       , pch = "."
       , cex = 4
       , log = "xy"
       , las = 1
       , xaxp = c(1,10,1)
       , yaxp = c(1,10,1)
       , labels = samples)
  
  # Return a correlation matrix
  corr.m <- matrix(1, nr.samples, nr.samples)
  colnames(corr.m) <- samples
  rownames(corr.m) <- samples
  corr.m[lower.tri(corr.m)] <- corr.v
  corr.m[upper.tri(corr.m)] <- t(corr.m)[upper.tri(corr.m)]
  corr.m
}

# Vanja's version of smooth scatter that allows passing range.x argument to grDevices:::.smoothScatterCalcDensity function to calculate 2D kernel smoothed density

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
