setGeneric(
name="plotCorrelation",
def=function(object, what = "CTSS", values = "raw", samples = "all", method = "pearson", tagCountThreshold = 1, applyThresholdBoth = FALSE, plotSize=800){
	standardGeneric("plotCorrelation")
}
)

setMethod("plotCorrelation",
signature(object = "CAGEset"),
function (object, what = "CTSS", values = "raw", samples = "all", method = "pearson", tagCountThreshold = 1, applyThresholdBoth = FALSE, plotSize=800){
	
	sample.labels <- sampleLabels(object)
	
	if(what == "CTSS"){
		if(values == "raw"){
			tag.count <- object@tagCountMatrix
		}else if(values == "normalized"){
			tag.count <- object@normalizedTpmMatrix
		}else{
			stop("'values' parameter must be one of the (\"raw\", \"normalized\")")
		}
	}else if(what == "consensusClusters"){
		tag.count <- as.data.frame(object@consensusClustersTpmMatrix)
	}else{
		stop("'what' parameter must be one of the (\"CTSS\", \"consensusClusters\")")
	}


	if(all(samples %in% sample.labels)){
		tag.count <- tag.count[,samples]
		nr.samples <- length(samples)
	}else if(samples == "all"){
		samples <- sample.labels
		nr.samples <- length(samples)
	}else{
		stop("'samples' parameter must be either \"all\" or a character vector of valid sample labels!")
	}
	
	corr.m <- matrix(rep(1, (nr.samples)^2), nrow = nr.samples)
	colnames(corr.m) <- samples
	rownames(corr.m) <- samples
	
	filename <- paste(what, "_", values, "_values_pairwise_correlation.png", sep = "")
	png(filename = filename, width = (plotSize + 36) * nr.samples + 0.1*plotSize*(log2(nr.samples)/log2(3) + (2-1/log2(3))), height = (plotSize + 36) * nr.samples + 0.05*plotSize*(log2(nr.samples)/log2(3) + (2-1/log2(3))), family = "Helvetica", res = 360)
	par(mfrow = c(nr.samples, nr.samples), mai = c(0.05,0.05,0.05,0.05), omi = c(0.05*plotSize*(log2(nr.samples)/log2(3) + (2-1/log2(3)))/360,0.1*plotSize*(log2(nr.samples)/log2(3) + (2-1/log2(3)))/360,0,0))
	
	for(i in c(1:nr.samples)){
		for(j in c(1:nr.samples)){
			
			if(i == j){
				
				plot(1, 1, type = "n", bty = "n", xlim = c(0,1), ylim = c(0,1), axes = F)
				text(0.5, 0.5, samples[i], cex = 0.98/max(sapply(samples, strwidth)))
				box(lwd = 3)
		
			}else if(j > i){
			
				x <- tag.count[,samples[j]]
				y <- tag.count[,samples[i]]
				if(applyThresholdBoth){
					idx <- (x >= tagCountThreshold) & (y >= tagCountThreshold)
				}else{
					idx <- (x >= tagCountThreshold) | (y >= tagCountThreshold)
				}
				x <- x[idx]
				y <- y[idx]
			
				pairwise.cor <- cor(x = x, y = y, method = method)
				plot(1, 1, type = "n", bty = "n", xlim = c(0,1), ylim = c(0,1), axes = F)
				txt <- sprintf("%.2f", pairwise.cor)
				txt.abs <- sprintf("%.2f", abs(pairwise.cor))
				text(0.5, 0.5, txt, cex = 1.5 + 0.5/strwidth(txt.abs) * abs(pairwise.cor))
				box(lwd = 3)
				corr.m[i,j] <- pairwise.cor
				corr.m[j,i] <- pairwise.cor
			
			}else{
			
				x <- tag.count[,samples[j]]
				y <- tag.count[,samples[i]]
				if(applyThresholdBoth){
					idx <- (x >= tagCountThreshold) & (y >= tagCountThreshold)
				}else{
					idx <- (x >= tagCountThreshold) | (y >= tagCountThreshold)
				}
				x <- x[idx]
				y <- y[idx]
				
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
				
#				fit.lm <- lm(y ~ x)
#				.mySmoothScatter(x = x, y = y, xlim = c(0, 200), ylim = c(0,200), nrpoints = 0, nbin = c(800, 800), bandwidth = c(3, 3), transformation = function(x) x^(1/9), axes = F)
#				lines(x = c(0,10,100,1000), y = coefficients(fit.lm)[2]*c(0,10,100,1000) + coefficients(fit.lm)[1], col = "red3", lwd = 3)				
				
			}
		}
	}
	
	dev.off()
	message("\nFile '", filename, "' has been created in your working directory (", getwd(), ")")
	
	return(corr.m)
	
}
)


# my version of smooth scatter that allows passing range.x argument to grDevices:::.smoothScatterCalcDensity function to calculate 2D kernel smoothed density

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
    map <- grDevices:::.smoothScatterCalcDensity(x, nbin, bandwidth, range.x = list(xlim = c(xlim[1] - 1.5*bandwidth[1], xlim[2] + 1.5*bandwidth[1]), ylim = c(ylim[1] - 1.5*bandwidth[2], ylim[2] + 1.5*bandwidth[2])))
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


