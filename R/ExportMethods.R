
#' @include CAGEr.R

###########################################################
# Functions for exporting results (graphical and textual)
#

#' @name plotReverseCumulatives
#' 
#' @title Plotting reverse cumulative number of CAGE tags per CTSS
#' 
#' @description Creates PDF file with plots of reverse cumulative number of CAGE
#' tags per CTSS for all CAGE datasets present in the \code{\link{CAGEr}} object.
#' The plots should be used as help in choosing the parameters for power-law
#' normalization: the range of values to fit the power-law and the slope of the
#' referent power-law distribution (Balwierz \emph{et al}., Genome Biology 2009).
#' 
#' @param object A CAGEr object
#' 
#' @param values Specifies which values should be plotted.  Can be either
#' \code{"raw"} to plot reverse cumulatives of raw CAGE tag counts or
#' \code{"normalized"} to plot normalized tag count values.
#' 
#' @param fitInRange An integer vector with two values specifying a range of tag
#' count values to be used for fitting a power-law distribution to reverse
#' cumulatives.  Used only when \code{values = "raw"}, otherwise ignored.
#' See Details.
#' 
#' @param onePlot Logical, should all CAGE datasets be plotted in the same
#' plot (TRUE) or in separate plots (FALSE) within the same PDF file. 
#' 
#' @details Number of CAGE tags (X-axis) is plotted against the number of TSSs that
#' are supported by >= of that number of tags (Y-axis) on a log-log scale for each
#' sample. In addition, a power-law distribution is fitted to each reverse cumulative
#' using the values in the range specified by \code{fitInRange} parameter. The fitted
#' distribution is defined by \code{y = -1 * alpha * x + beta} on the log-log scale,
#' and the value of \code{alpha} for each sample is shown on the plot. In addition,
#' a suggested referent power-law distribution to which all samples should be
#' normalized is drawn on the plot and corresponding parameters (slope alpha and total
#' number of tags T) are denoted on the plot.  Referent distribution is chosen so
#' that its slope (alpha) is the median of slopes fitted to individual samples and
#' its total number of tags (T) is the power of 10 nearest to the median number of tags
#' of individual samples.  Resulting plots are helpful in deciding whether power-law
#' normalization is appropriate for given samples and reported \code{alpha} values aid
#' in choosing optimal \code{alpha} value for referent power-law distribution to
#' which all samples will be normalized. For details about normalization see
#' \code{\link{normalizeTagCount}} function.
#' 
#' @return Plots of reverse cumulative number of CAGE tags per CTSS for each CAGE
#' dataset within CAGEr object.  Alpha values of fitted power-laws and suggested
#' referent power-law distribution are reported on the plot in case \code{values = "raw"}.
#' 
#' @references Balwierz \emph{et al}. (2009) Methods for analyzing deep sequencing
#' expression data: constructing the human and mouse promoterome with deepCAGE data,
#' \emph{Genome Biology} \bold{10}(7):R79. 
#' 
#' @author Vanja Haberle
#' 
#' @seealso \code{\link{normalizeTagCount}}
#' @family CAGEr plot functions
#' 
#' @examples 
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' plotReverseCumulatives(exampleCAGEset, values = "raw",fitInRange = c(10,500), onePlot = TRUE)
#' 
#' @importFrom VGAM zeta
#' @export

setGeneric(
name="plotReverseCumulatives",
def=function(object, values = "raw", fitInRange = c(10, 1000), onePlot = FALSE){
	standardGeneric("plotReverseCumulatives")
}
)

setMethod("plotReverseCumulatives",
signature(object = "CAGEr"),
function (object, values = "raw", fitInRange = c(10, 1000), onePlot = FALSE){
		
	sample.labels <- sampleLabels(object)
	if(values == "raw"){
		tag.count <- CTSStagCountDf(object)
	}else if(values == "normalized"){
		tag.count <- CTSSnormalizedTpmDf(object)
	}else{
		stop("'values' parameter must be one of the (\"raw\", \"normalized\")")
	}
	
	#pdf(file = paste("CTSS_reverse_cumulatives_", values, "_all_samples.pdf", sep = ""), width = 8, height = 8, onefile = T, bg = "transparent", family = "Helvetica", fonts = NULL)
	par(mar = c(5,5,5,2))
	cols <- names(sample.labels)
	
	if(values == "raw"){
		fit.coefs.m <- apply(tag.count, 2, function(x) {.fit.power.law.to.reverse.cumulative(values = as.integer(x), val.range = fitInRange)})
		fit.slopes <- fit.coefs.m[1,]
		names(fit.slopes) <- sample.labels
		reference.slope <- min(median(fit.slopes), -1.05)
		library.sizes <- librarySizes(object)
		reference.library.size <- 10^floor(log10(median(library.sizes)))
#reference.intercept <- log(reference.library.size/zeta(-1*reference.slope))  # intercept on natural logarithm scale
		reference.intercept <- log10(reference.library.size/VGAM::zeta(-1*reference.slope))  # intercept on log10 scale used for plotting with abline
	}else if(values == "normalized"){
#		fit.coefs.m <- apply(tag.count, 2, function(x) {.fit.power.law.to.reverse.cumulative(values = x, val.range = fitInRange)})
	}
	
	if(onePlot == TRUE){
		vals <- tag.count[, 1]		
		if(values == "raw"){
			.plotReverseCumulative(values = as.integer(vals), col = cols[1], title = "All samples")
			if(length(sample.labels) > 1){
				sapply(c(2:length(sample.labels)), function(x) {vals <- as.integer(tag.count[, sample.labels[x]]); .plotReverseCumulative(values = vals, col = cols[x], add = TRUE)})
			}
			abline(v = fitInRange, lty = "dotted")
			abline(a = reference.intercept, b = reference.slope, col = "#7F7F7F7F", lty = "longdash")
			legend("topright", legend = paste("(", formatC(-1*fit.slopes, format = "f", digits = 2), ") ", sample.labels, sep = ""), bty = "n", col = cols, text.col = cols, lwd = 2, cex = 1.3, y.intersp = 1.2)
			legend("bottomleft", legend = c("Ref. distribution:", paste(" alpha = ", sprintf("%.2f", -1*reference.slope), sep = ""), paste(" T = ", reference.library.size, sep = "")), bty = "n", col = NA, text.col = "#7F7F7F", cex = 1.3, y.intersp = 1.2)
		}else if(values == "normalized"){
			.plotReverseCumulative(values = vals, col = cols[1], title = "All samples")
			if(length(sample.labels) > 1){			
				sapply(c(2:length(sample.labels)), function(x) {vals <- tag.count[, sample.labels[x]]; .plotReverseCumulative(values = vals, col = cols[x], add = TRUE)})
			}
			legend("topright", legend = sample.labels, bty = "n", col = cols, text.col = cols, lwd = 2, cex = 1.3, y.intersp = 1.2)
		}
	}else{
		if(values == "raw"){
			sapply(sample.labels, function(x) {vals <- as.integer(tag.count[, x]); .plotReverseCumulative(values = vals, col = cols[which(sample.labels == x)], title = x, col.title = cols[which(sample.labels == x)]); abline(v = fitInRange, lty = "dotted"); abline(a = reference.intercept, b = reference.slope, col = "#7F7F7F7F", lty = "longdash"); text(min(fitInRange), 10^6, labels = paste(" alpha =", formatC(-1*fit.slopes[x], format = "f", digits = 2), sep = " "), adj = c(0,1), col = cols[which(sample.labels == x)], cex = 1.3); legend("bottomleft", legend = c("Ref. distribution:", paste(" alpha = ", sprintf("%.2f", -1*reference.slope), sep = ""), paste(" T = ", reference.library.size, sep = "")), bty = "n", col = NA, text.col = "#7F7F7F", cex = 1.3, y.intersp = 1.2)})
		}else if(values == "normalized"){
			sapply(sample.labels, function(x) {vals <- tag.count[, x]; .plotReverseCumulative(values = vals, col = cols[which(sample.labels == x)], title = x, col.title = cols[which(sample.labels == x)])})
		}
	}
	#dev.off()
	#message("\nFile 'CTSS_reverse_cumulatives_", values, "_all_samples.pdf' has been created in your working directory (", getwd(), ")")
	
}
)

#' @name exportCTSStoBedGraph
#' 
#' @title Creating bedGraph/bigWig tracks of CAGE transcription starts sites
#' 
#' @description Creates bedGraph or BigWig file(s) with track(s) of CAGE signal supporting
#' each TSS that can be visualised in the UCSC Genome Browser.
#' 
#' @param object A \code{\link{CAGEr}} object.
#' 
#' @param values Specifies which values will be exported to the bedGraph file. Can be either
#'        \code{"raw"} to export raw tag count values or \code{"normalized"} to export
#'        normalized values.
#' 
#' @param format The format of the output. 
#' 
#' @param oneFile Logical, should all CAGE datasets be exported as individual tracks into the
#'        same bedGraph file (TRUE) or into separate bedGraph files (FALSE). Used only when
#'        \code{format="bedGraph"}, otherwise ignored.
#' 
#' @return Creates bedGraph or BigWig file(s) in the working directory that can be directly
#' visualised as custom tracks in the UCSC Genome Browser.  If \code{format="bedGraph"} and
#' \code{oneFile = TRUE} one bedGraph file containing multiple annotated tracks will be created,
#' otherwise two files per CAGE dataset will be created, one for plus strand and one for minus
#' strand CTSSs, and they will be named according to the labels of individual datasets.  All
#' bedGraph files contain headers with track description and can be directly uploaded as custom
#' tracks to the UCSC Genome Browser. 
#' 
#' When \code{format="bigWig"}, two binary BigWig files per CAGE dataset are created, one for
#' plus strand and one for minus strand CTSSs. Since BigWig files cannot contain headers with
#' track description, a separate file named "CTSS.normalized.all.samples.track.description.txt"
#' is created, which contains track headers for all BigWig files. To use these headers for
#' adding custom tracks to the UCSC Genome Browser, move the BigWig files to a web location and
#' edit the bigDataUrl sections in the headers file to point to corresponding BigWig files.
#' 
#' @author Vanja Haberle
#' 
#' @family CAGEr export functions
#' 
#' @seealso \code{\link{normalizeTagCount}}
#' 
#' @examples
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' exportCTSStoBedGraph(exampleCAGEset, values = "normalized", format = "bedGraph", oneFile = TRUE)
#' 
#' @export

setGeneric(
name="exportCTSStoBedGraph",
def=function(object, values = "normalized", format = "BigWig", oneFile = TRUE){
	standardGeneric("exportCTSStoBedGraph")
}
)

setMethod("exportCTSStoBedGraph",
signature(object = "CAGEr"),
function (object, values = "normalized", format = "BigWig", oneFile = TRUE){

  if (values == "raw") {
    data <- SummarizedExperiment( rowRanges = CTSScoordinatesGR(object)
                                , assay     = SimpleList(CTSStagCountDF(object)))
  } else if (values == "normalized"){
    data <- SummarizedExperiment( rowRanges = CTSScoordinatesGR(object)
                                , assay     = SimpleList(CTSSnormalizedTpmDF(object)))
  } else {
    stop("'values' parameter must be one of the (\"raw\", \"normalized\")")
  }

  if (format == "BigWig"){
    genome <- getRefGenome(genomeName(object))
    .export.bw.all(data = data, sample.labels = sampleLabels(object), v = values, genome = genome)
  } else if (format == "bedGraph"){
    .export.bedgraph.all(data = data, sample.labels = sampleLabels(object), v = values, oneFile = oneFile)
  } else {
    stop("'format' parameter must be one of the (\"BigWig\", \"bedGraph\")")
  }
	
  message("\n", format, " file(s) for CTSS ", values, " counts have been created in your working directory (", getwd(), ")")
  invisible(1)
})

#' @name plotInterquantileWidth
#' @noRd
#' @export

setGeneric(
name="plotInterquantileWidth",
def=function(object, clusters, tpmThreshold = 5, qLow = 0.1, qUp = 0.9, xlim = c(0,150), ...){
	standardGeneric("plotInterquantileWidth")
}
)

setMethod("plotInterquantileWidth",
signature(object = "CAGEset"),
function (object, clusters, tpmThreshold, qLow, qUp, xlim = c(0,150), ...){
	
	sample.labels <- sampleLabels(object)
	cols <- names(sample.labels)
	
	if(clusters == "tagClusters"){	
		
		if(length(object@tagClustersQuantileLow)>0 & length(object@tagClustersQuantileUp)>0) {
			if(!(paste("q_", qLow, sep = "") %in% colnames(object@tagClustersQuantileLow[[1]]) & paste("q_", qUp, sep = "") %in% colnames(object@tagClustersQuantileUp[[1]]))){
				stop("No data for given quantile positions! Run 'quantilePositions()' function for desired quantiles first!")
			}
		}else{
			stop("No data for given quantile positions! Run 'quantilePositions()' function for desired quantiles first!")
		}
		
		filename <- "TC"
		q.low <- object@tagClustersQuantileLow
		q.up <- object@tagClustersQuantileUp
		idx.list <- lapply(as.list(sample.labels), function(x) {
			   
								cl <- tagClusters(object, sample = x)
								idx <- cl$tpm >= tpmThreshold
								return(idx)
						   
							}
						   )
	
	}else if (clusters == "consensusClusters"){
		
		if(length(object@consensusClustersQuantileLow)>0 & length(object@consensusClustersQuantileUp)>0) {
			
			if(!(paste("q_", qLow, sep = "") %in% colnames(object@consensusClustersQuantileLow[[1]]) & paste("q_", qUp, sep = "") %in% colnames(object@consensusClustersQuantileUp[[1]]))){
				stop("No data for given quantile positions! Run 'quantilePositions()' function for desired quantiles first!")
			}
		}else{
			stop("No data for given quantile positions! Run 'quantilePositions()' function for desired quantiles first!")
		}
		
		filename <- "consensusClusters_interquantile_width_all_samples.pdf"
		q.low <- object@consensusClustersQuantileLow
		q.up <- object@consensusClustersQuantileUp
		cl <- object@consensusClustersTpmMatrix
		idx.list <- lapply(as.list(sample.labels), function(x) {idx <- cl[,x][cl[,x] > 0] >= tpmThreshold})		
		
	}else{
		stop("'clusters' parameter must be one of the (\"tagClusters\", \"consensusClusters\")")
	}

	names(idx.list) <- sample.labels
	
	pdf(file = paste(clusters, "_interquantile_width_all_samples.pdf", sep = ""), width = 8, height = 8, onefile = T, bg = "transparent", family = "Helvetica", fonts = NULL)
	par(mar = c(5,5,5,1))
	sapply(sample.labels, function(x) {
		   
		   q.low.s <- q.low[[x]]
		   q.low.s <- as.integer(q.low.s[, which(colnames(q.low.s) == paste("q_", qLow, sep = ""))])
		   q.up.s <- q.up[[x]]
		   q.up.s <- as.integer(q.up.s[, which(colnames(q.up.s) == paste("q_", qUp, sep = ""))])
		   width <- q.up.s[idx.list[[x]]] - q.low.s[idx.list[[x]]] + 1
		   h <- hist(width, breaks = round(max(width)/2), plot = F)
		   h$counts <- h$counts/sum(h$counts)
		   col <- as.integer(col2rgb(cols[which(sample.labels == x)]))/255
		   plot(h, xlim = xlim, main = x, xlab = paste(clusters, " interquantile width q", qLow, "-q", qUp, " (bp)", sep = ""), ylab = "relative frequency", col = rgb(col[1], col[2], col[3], 0.5), border = cols[which(sample.labels == x)], cex.axis = 1.8, cex.lab = 1.8, cex.main = 2.5, col.main = cols[which(sample.labels == x)], ...)
		   
		   }
		   )
	dev.off()
	message("\nFile '", clusters, "_interquantile_width_all_samples.pdf' has been created in your working directory (", getwd(), ")")
	
}
)

#' @name plotExpressionProfiles
#' @noRd
#' @export

setGeneric(
name="plotExpressionProfiles",
def=function(object, what){
	standardGeneric("plotExpressionProfiles")
}
)

setMethod("plotExpressionProfiles",
signature(object = "CAGEset"),
function (object, what){
	
	if(what == "CTSS") {
		
		cl <- object@CTSSexpressionClasses
		if(length(cl)>0){
			tpm.mx <- object@normalizedTpmMatrix
			tpm.mx <- tpm.mx[as.integer(names(cl)),]
			cl.method <- object@CTSSexpressionClusteringMethod
		}else{
			stop("No CTSS expression profiling has been done yet! Run getExpressionProfiles function first!")
		}
		
	}else if(what == "consensusClusters"){
		
		cl <- object@consensusClustersExpressionClasses
		if(length(cl)>0){
			tpm.mx <- object@consensusClustersTpmMatrix
			tpm.mx <- tpm.mx[as.integer(names(cl)),]
			cl.method <- object@consensusClustersExpressionClusteringMethod
		}else{
			stop("No consensusClusters expression profiling has been done yet! Run getExpressionProfiles function first!")
		}
		
	}else{
		stop("'what' parameter must be one of the (\"CTSS\", \"consensusClusters\")")
	}
	
	l <- .extract.cluster.info(tpm.mx, cl)
	cl <- l[[1]]
	m <- l[[2]]
	file_name = paste(what, "_expression_profiles.pdf", sep = "")
	pdf(file = file_name, height = 5.5 * (max(cl[,2]) + 1) + 3, width = 6 * (max(cl[,1]) + 1))
	par(omi = c(3,0.5,0.5,0.5), mfrow = c(max(cl[,2]) + 1, max(cl[,1]) + 1))
	suppressWarnings(.plot.clusters.beanplots(value.matrix = m, cl = cl, cl.method = cl.method, dim.som.x = max(cl[,1]) + 1, dim.som.y = max(cl[,2]) + 1, cex.axis = 3, las = 2))
	dev.off()
	message("\nFile '", file_name, "' has been created in your working directory (", getwd(), ")")	
	
}
)

#' @name exportToBed
#' @noRd
#' @export

setGeneric(
name="exportToBed",
def=function(object, what, qLow = NULL, qUp = NULL, colorByExpressionProfile = FALSE, oneFile = TRUE){
	standardGeneric("exportToBed")
}
)


setMethod("exportToBed",
signature(object = "CAGEset"),
function (object, what, qLow = NULL, qUp = NULL, colorByExpressionProfile = FALSE, oneFile = TRUE){
	
	sample.labels <- sampleLabels(object)

	if(what == "CTSS") {
		
		oneFile <- TRUE
		use.blocks <- F
		ctss <- object@CTSScoordinates
		#filtered_ctss <- object@filteredCTSSidx

		if(colorByExpressionProfile == TRUE){
			cl <- object@CTSSexpressionClasses
			n <- names(cl)
			cl <- .extract.cluster.info(cl = cl)
			cl <- data.frame(ctss = n, x.cor = cl[,1], y.cor = cl[,2])		
			ctss <- merge(cl, ctss, by.x = "ctss", by.y = 0, all.x = T, all.y = F)
			ctss <- data.frame(ctss = ctss$ctss, chr = ctss$chr, start = ctss$pos-1, end = ctss$pos, strand = ctss$strand, x.cor = ctss$x.cor, y.cor = ctss$y.cor)
			track.file <- "CTSS.colored.by.expression.profile.bed"
			track.names <- list("CTSS (colored by expression profile)")
			clustering.method <- object@CTSSexpressionClusteringMethod
		}else{
			ctss <- data.frame(chr = ctss$chr, start = ctss$pos-1, end = ctss$pos, strand = ctss$strand)
			track.file <- "CTSS.pooled.samples.bed"
			track.names <- list("CTSS (pooled samples)")
			#filtered_cols <- c("TRUE" = c("0,0,0"), "FALSE" = c("127,127,127"))
			#cols = list(filtered_cols[as.character(filtered_ctss)])
			cols = list(c("0,0,0"))
		}
		clusters.q.list = list(ctss)
		
	}else if(what == "tagClusters") {
		
		colorByExpressionProfile <- FALSE
		
		if(length(qLow) > 0 & (length(object@tagClustersQuantileLow)>0 & length(object@tagClustersQuantileLow)>0)) {
		
		if(paste("q_", qLow, sep = "") %in% colnames(object@tagClustersQuantileLow[[1]]) & paste("q_", qUp, sep = "") %in% colnames(object@tagClustersQuantileUp[[1]])){
		
		use.blocks <- T
		q.low <- object@tagClustersQuantileLow
		q.low <- lapply(q.low, function(x) {colnames(x)[2:ncol(x)] <- paste("qLow_", do.call(rbind, strsplit(colnames(x)[2:ncol(x)], split = "_", fixed = T))[,2], sep = ""); return(x)})
		q.up <- object@tagClustersQuantileUp
		q.up <- lapply(q.up, function(x) {colnames(x)[2:ncol(x)] <- paste("qUp_", do.call(rbind, strsplit(colnames(x)[2:ncol(x)], split = "_", fixed = T))[,2], sep = ""); return(x)})
		q <- lapply(as.list(1:length(sample.labels)), function(x) {merge(q.low[[x]], q.up[[x]], by.x = "cluster", by.y = "cluster")})
		names(q) <- sample.labels
		clusters <- lapply(as.list(sample.labels), function(x) {tagClusters(object, sample = x)})
		names(clusters) <- sample.labels
		clusters.q.list <- lapply(as.list(sample.labels), function(x) {merge(clusters[[x]], q[[x]], by.x = "cluster", by.y = "cluster")})
		track.names <- paste(sample.labels, paste(" (tag clusters (TC) q(", qLow, ")-q(",qUp,"))", sep = ""), sep = "")
		r <- paste(".qLow", qLow, "_qUp", qUp, sep = "")
			
		}else{
			stop("No data for given quantile positions! Run 'quantilePositions()' function for desired quantiles first, or omit 'qLow' and 'qUp' parameters to use start and end coordinates instead!")
		}
		}else if(length(qLow) > 0){
			stop("No data for given quantile positions! Run 'quantilePositions()' function for desired quantiles first, or omit 'qLow' and 'qUp' parameters to use start and end coordinates instead!")
		}else{
			use.blocks <- F
			clusters.q.list <- lapply(as.list(sample.labels), function(x) {tagClusters(object, sample = x)})			
			track.names <- paste(sample.labels, paste(" (tag clusters (TC))", sep = ""), sep = "")
			r <- ""
		}
		
		names(clusters.q.list) <- sample.labels
		itemRgb = FALSE
		cols <- names(sample.labels)
		cols <- as.list(apply(sapply(cols, function(x) {as.integer(col2rgb(x))}), 2, function(y) {paste(y, collapse = ",")}))
		names(cols) <- sample.labels

		if(oneFile){
			track.file <- rep(paste("All.samples.tagClusters", r, ".bed", sep = ""), length(clusters.q.list))
		}else{
			track.file <- paste(sample.labels, ".tagClusters", r, ".bed", sep = "")
		}
		
	}else if(what == "consensusClusters"){
		
		clusters <- object@consensusClusters
		colnames(clusters)[1] = "cluster"
		if(!(colorByExpressionProfile)){
			cols <- as.list(rep("0,0,0", length(sample.labels)))
		}else{
			clustering.method <- object@consensusClustersExpressionClusteringMethod		
		cl <- object@consensusClustersExpressionClasses
		n <- names(cl)
		cl <- .extract.cluster.info(cl = cl)
		cl <- data.frame(cluster = n, x.cor = cl[,1], y.cor = cl[,2])		
		clusters <- merge(clusters, cl, by.x = "cluster", by.y = "cluster")
		}
		
		if(length(qLow) > 0 & (length(object@consensusClustersQuantileLow)>0 & length(object@consensusClustersQuantileUp)>0)) {

		if(paste("q_", qLow, sep = "") %in% colnames(object@consensusClustersQuantileLow[[1]]) & paste("q_", qUp, sep = "") %in% colnames(object@consensusClustersQuantileUp[[1]])){
		
		use.blocks <- T
		q.low <- object@consensusClustersQuantileLow
		q.low <- lapply(q.low, function(x) {colnames(x)[2:ncol(x)] <- paste("qLow_", do.call(rbind, strsplit(colnames(x)[2:ncol(x)], split = "_", fixed = T))[,2], sep = ""); return(x)})
		q.up <- object@consensusClustersQuantileUp
		q.up <- lapply(q.up, function(x) {colnames(x)[2:ncol(x)] <- paste("qUp_", do.call(rbind, strsplit(colnames(x)[2:ncol(x)], split = "_", fixed = T))[,2], sep = ""); return(x)})
		q <- lapply(as.list(1:length(sample.labels)), function(x) {merge(q.low[[x]], q.up[[x]], by.x = "cluster", by.y = "cluster")})
		names(q) <- sample.labels
						
		cumsums <- object@CTSScumulativesConsensusClusters
		dom.pos <- list()
		for(i in 1:length(cumsums)){
			a <- lapply(cumsums[[i]], function(y) {.get.dominant.ctss(as.numeric(y), isCumulative = T)})
			b <- data.frame(cluster = as.integer(names(a)), dominant_ctss = unlist(a))
			dom.pos[[i]] <- b
		}
		names(dom.pos) <- names(cumsums)
		
		clusters.q.list <- lapply(as.list(sample.labels), function(x) {a <- merge(clusters, dom.pos[[x]], by.x = "cluster", by.y = "cluster", all.x = F, all.y = T); a$dominant_ctss <- a$start + a$dominant_ctss; b <- merge(a, q[[x]], by.x = "cluster", by.y = "cluster"); return(b)})
		track.names <- paste(sample.labels, paste(" (consensus clusters q(", qLow, ")-q(",qUp,"))", sep = ""), sep = "")
		r <- paste(".qLow", qLow, "_qUp", qUp, sep = "")
		quantiles <- T
		}else{
			stop("No data for given quantile positions! Run 'quantilePositions()' function for desired quantiles first, or omit 'qLow' and 'qUp' parameters to use start and end coordinates instead!")
		}
		}else if(length(qLow) > 0){
			stop("No data for given quantile positions! Run 'quantilePositions()' function for desired quantiles first, or omit 'qLow' and 'qUp' parameters to use start and end coordinates instead!")
		}else{
			use.blocks <- F
			clusters.q.list <- list(clusters)			
			track.names <- "Consensus clusters"
			r <- ""
			oneFile <- T
			quantiles <- F
		}
				   
		if(oneFile){
			if(quantiles){
			track.file <- rep(paste("All.samples.consensusClusters", r, ".bed", sep = ""), length(clusters.q.list))
			}else{
				track.file <- "ConsensusClusters.bed"
			}
		}else{
			track.file <- paste(sample.labels, ".consensusClusters.", r, ".bed", sep = "")		
		}
		
	}else{
		stop("'what' parameter must be one of the (\"CTSS\", \"tagClusters\", \"consensusClusters\")")
	}

	if(colorByExpressionProfile == TRUE){
		
		itemRgb <- TRUE
		cols.init <- c("red", "gold", "green", "blue")
		color.matrix.solid <- .myColorMatrix(matrix(cols.init, nrow = 2), nrows = max(cl$y.cor)+1, ncols = max(cl$x.cor)+1)
		color.matrix.solid <- matrix(as.vector(color.matrix.solid), ncol = max(cl$x.cor)+1, byrow = F)
		cols <- lapply(clusters.q.list, function(x) {m <- as.matrix(x[,c("y.cor", "x.cor")]); a <- apply(t(col2rgb(color.matrix.solid[m+1])), 1, function(x) {paste(x, collapse = ',')}); return(a)})
		if(clustering.method == "som"){
			names <- lapply(clusters.q.list, function(x) {paste(x[,"x.cor"], x[,"y.cor"], sep = "_")})
		}else if(clustering.method == "kmeans"){
			names <- lapply(clusters.q.list, function(x) {x[,"x.cor"]})
		}
		
	}else{
		itemRgb <- FALSE
		names <- as.list(rep(".", length(clusters.q.list)))
	}
	
	track.descriptions <- track.names
	
	for(i in 1:length(clusters.q.list)){
		if(i == 1){
			app <- F
		}else{
			app <- T
		}
		if(!oneFile){
			if(file.exists(track.file[i])){
				file.remove(track.file[i])
			}				
		}
		.make.cluster.bed.track(clusters.q = clusters.q.list[[i]], use.blocks = use.blocks, q.low = qLow, q.up = qUp, track.file = track.file[i], track.name = track.names[i], track.description = track.descriptions[i], cols = cols[[i]], name = names[[i]], itemRgb = itemRgb, app = app)
		
	}
	
	if(oneFile){
		message("\nFile '", track.file[1], "' has been created in your working directory (", getwd(), ")")	
	}else{
		message("\nFiles '", sub(sample.labels[1], "*", track.file[1]), "' for all samples have been created in your working directory (", getwd(), ")")
	}
	
}
)





