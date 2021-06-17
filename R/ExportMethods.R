#' @include CAGEr.R

###########################################################
# Functions for exporting results (graphical and textual)
#

#' @name plotReverseCumulatives
#' 
#' @title Plot reverse cumulative number of CAGE tags per CTSS
#' 
#' @description Plots the reverse cumulative distribution of the number of CAGE
#' tags per CTSS for all CAGE datasets present in the [`CAGEr`] object.
#' The plots should be used as help in choosing the parameters for power-law
#' normalization: the range of values to fit the power-law and the slope of the
#' referent power-law distribution (Balwierz _et al_., Genome Biology 2009).
#' 
#' @param object A CAGEr object
#' 
#' @param values Which values should be plotted: `raw` (default) for raw CAGE
#' tag counts or `normalized` for normalized tag count values.
#' 
#' @param fitInRange An integer vector with two values specifying a range of tag
#' count values to be used for fitting a power-law distribution to reverse
#' cumulatives.  Ignored is set to `NULL`.  See Details.
#' 
#' @param onePlot Logical, should all CAGE datasets be plotted in the same
#' plot (`TRUE`) or in separate plots (`FALSE`). 
#' 
#' @param main Main title for the plot.
#' 
#' @param legend Set to `NULL` to prevent the display of the sample legend.
#' 
#' @param xlab,ylab Axis labels passed to [`plot`].
#' 
#' @param xlim,ylim Axis range parameters passed to [`plot`].
#' 
#' @details Number of CAGE tags (X-axis) is plotted against the number of TSSs that
#' are supported by >= of that number of tags (Y-axis) on a log-log scale for each
#' sample. In addition, a power-law distribution is fitted to each reverse cumulative
#' using the values in the range specified by `fitInRange` parameter. The fitted
#' distribution is defined by `y = -1 * alpha * x + beta` on the log-log scale,
#' and the value of `alpha` for each sample is shown on the plot. In addition,
#' a suggested referent power-law distribution to which all samples should be
#' normalized is drawn on the plot and corresponding parameters (slope alpha and total
#' number of tags T) are denoted on the plot.  Referent distribution is chosen so
#' that its slope (alpha) is the median of slopes fitted to individual samples and
#' its total number of tags (T) is the power of 10 nearest to the median number of tags
#' of individual samples.  Resulting plots are helpful in deciding whether power-law
#' normalization is appropriate for given samples and reported `alpha` values aid
#' in choosing optimal `alpha` value for referent power-law distribution to
#' which all samples will be normalized. For details about normalization see
#' [`normalizeTagCount`] function.
#' 
#' @return Plots of reverse cumulative number of CAGE tags per CTSS for each CAGE
#' dataset within CAGEr object.  Alpha values of fitted power-laws and suggested
#' referent power-law distribution are reported on the plot in case `values = "raw"`.
#' 
#' @references Balwierz _et al_. (2009) Methods for analyzing deep sequencing
#' expression data: constructing the human and mouse promoterome with deepCAGE data,
#' _Genome Biology_ **10**(7):R79. 
#' 
#' @author Vanja Haberle
#' 
#' @seealso [`normalizeTagCount`]
#' @family CAGEr plot functions
#' 
#' @examples 
#' 
#' plotReverseCumulatives( exampleCAGEexp, xlim = c(1, 1e4), ylim = c(1, 1e5)
#'                       , fitInRange = c(5,100), onePlot = TRUE)
#' plotReverseCumulatives( exampleCAGEexp, values = "normalized"
#'                       , fitInRange = c(200, 2000), onePlot = TRUE)
#' plotReverseCumulatives( exampleCAGEexp[,4:5], fitInRange = c(5,100)
#'                       , onePlot = TRUE, main = "prim6 replicates")
#' 
#' @importFrom graphics abline box legend par plot strwidth text
#' @importFrom stats cor median
#' @importFrom VGAM zeta
#' @export

setGeneric( "plotReverseCumulatives"
          , function( object, values = c("raw", "normalized")
                    , fitInRange = c(10, 1000)
                    , onePlot = FALSE, main = NULL, legend = TRUE
                    , xlab = "number of CAGE tags"
                    , ylab = "number of CTSSs (>= nr tags)"
                    , xlim = c(1, 1e5)
                    , ylim = c(1, 1e6))
	standardGeneric("plotReverseCumulatives"))

#' @rdname plotReverseCumulatives

setMethod( "plotReverseCumulatives", "CAGEr"
         , function ( object, values, fitInRange, onePlot, main, legend
                    , xlab, ylab, xlim, ylim) {
	sample.labels <- sampleLabels(object)
	values <- match.arg(values)
	#pdf(file = paste("CTSS_reverse_cumulatives_", values, "_all_samples.pdf", sep = ""), width = 8, height = 8, onefile = T, bg = "transparent", family = "Helvetica", fonts = NULL)
	old.par <- par(mar = c(5,5,5,2))
	on.exit(par(old.par))
	cols <- names(sample.labels)
	
  tag.count <- switch( values
                     , raw        = CTSStagCountDF(object)
                     , normalized = CTSSnormalizedTpmDF(object))

	if(! is.null(fitInRange)) {
		fit.coefs.m <- as.matrix(data.frame(lapply(tag.count, function(x) {.fit.power.law.to.reverse.cumulative(values = decode(x), val.range = fitInRange)})))
		fit.slopes <- fit.coefs.m[1,]
		names(fit.slopes) <- sample.labels
		reference.slope <- min(median(fit.slopes), -1.05)
		reference.library.size <- 10^floor(log10(median(sapply(tag.count, sum))))
#reference.intercept <- log(reference.library.size/zeta(-1*reference.slope))  # intercept on natural logarithm scale
		reference.intercept <- log10(reference.library.size/VGAM::zeta(-1*reference.slope))  # intercept on log10 scale used for plotting with abline
	}
	
	if(onePlot == TRUE){
		.plotReverseCumulative( values = tag.count[, 1]
		                      , col    = cols[1]
		                      , title  = ifelse(is.null(main), "All samples", main)
		                      , xlab   = xlab, ylab = ylab
		                      , xlim   = xlim, ylim = ylim)
		if(length(sample.labels) > 1){
			sapply(c(2:length(sample.labels)), function(x)
			  .plotReverseCumulative( values = tag.count[, sample.labels[x]]
			                        , col = cols[x]
			                        , add = TRUE))
		}
		if(!is.null(fitInRange)) {
			abline(v = fitInRange, lty = "dotted")
			abline(a = reference.intercept, b = reference.slope, col = "#7F7F7F7F", lty = "longdash")
			main.legend.text <- sprintf("(%.2f) %s", -1 * fit.slopes, sample.labels)
			legend("bottomleft", legend = c("Ref. distribution:", paste(" alpha = ", sprintf("%.2f", -1*reference.slope), sep = ""), paste(" T = ", reference.library.size, sep = "")), bty = "n", col = NA, text.col = "#7F7F7F", cex = 1.3, y.intersp = 1.2)
		} else {
			main.legend.text <- sample.labels
		}
		if (isTRUE(legend))
		  legend( "topright"
		        , legend    = main.legend.text
		        , bty       = "n"
		        , col       = cols
		        , text.col  = cols
		        , lwd       = 2
		        , cex       = 1.3
		        , y.intersp = 1.2)
	}else{
		if(!is.null(fitInRange)) {
			sapply(sample.labels, function(x) {
				.plotReverseCumulative( values    = tag.count[, x]
				                      , col       = cols[which(sample.labels == x)]
				                      , title     = x
				                      , col.title = cols[which(sample.labels == x)]
				                      , xlab      = xlab, ylab = ylab
				                      , xlim      = xlim, ylim = ylim)
				abline(v = fitInRange, lty = "dotted"); abline(a = reference.intercept, b = reference.slope, col = "#7F7F7F7F", lty = "longdash"); text(min(fitInRange), 10^6, labels = paste(" alpha =", formatC(-1*fit.slopes[x], format = "f", digits = 2), sep = " "), adj = c(0,1), col = cols[which(sample.labels == x)], cex = 1.3); legend("bottomleft", legend = c("Ref. distribution:", paste(" alpha = ", sprintf("%.2f", -1*reference.slope), sep = ""), paste(" T = ", reference.library.size, sep = "")), bty = "n", col = NA, text.col = "#7F7F7F", cex = 1.3, y.intersp = 1.2)})
		} else {
			sapply( sample.labels, function(x)
				.plotReverseCumulative( values    = tag.count[, x]
				                      , col       = cols[which(sample.labels == x)]
				                      , title     = x
				                      , col.title = cols[which(sample.labels == x)]
				                      , xlab      = xlab, ylab = ylab
				                      , xlim      = xlim, ylim = ylim)
			)
		}
	}
	invisible(TRUE)
})


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
#' exportCTSStoBedGraph(exampleCAGEexp, values = "normalized", format = "bedGraph", oneFile = TRUE)
#' exportCTSStoBedGraph(exampleCAGEexp, values = "normalized", format = "BigWig", oneFile = TRUE)
#' 
#' @export

setGeneric( "exportCTSStoBedGraph"
          , function(object, values = "normalized", format = "BigWig", oneFile = TRUE)
            	standardGeneric("exportCTSStoBedGraph"))

#' @rdname exportCTSStoBedGraph

setMethod("exportCTSStoBedGraph", "CAGEr", function (object, values, format, oneFile) {
  rr <- CTSScoordinatesGR(object)
  genome <- getRefGenome(genomeName(object))
  seqinfo(rr) <- seqinfo(genome)[seqlevels(rr)]

  if (values == "raw") {
    data <- SummarizedExperiment( rowRanges = rr
                                , assays    = SimpleList(CTSStagCountDF(object)))
  } else if (values == "normalized"){
    data <- SummarizedExperiment( rowRanges = rr
                                , assays    = SimpleList(CTSSnormalizedTpmDF(object)))
  } else {
    stop("'values' parameter must be one of the (\"raw\", \"normalized\")")
  }

  if (format == "BigWig"){
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
#' 
#' @title Plot cluster widths
#' 
#' @description Histograms of the interquantile width of tag clusters or
#' consensus clusters in each CAGE experiment.
#' 
#' @param object A [`CAGEexp`] object
#' 
#' @param clusters `tagClusters` or `consensusClusters`.
#' 
#' @param tpmThreshold Exclude clusters with normalized signal lower than `tpmThreshold`.
#' 
#' @param qLow,qUp Quantile defining the 5' ("lower") and 3' ("upper")
#'        boundaries of the clusters.
#' 
#' @param xlim Maximal width to be plotted.
#' 
#' @details Interquantile width is a more robust measure of the promoter width
#' than the total span of the region, because it takes into account the
#' magnitude of the expression in the region.  Positions of specified quantiles
#' within each cluster have to be calculated beforehand by calling
#' [`quantilePositions`].
#' 
#' @return Plots the histograms with the `ggplot2` engine and returns the plot
#' object invisibly.
#' 
#' @author Vanja Haberle
#' @author Charles Plessy
#' 
#' @family CAGEr plot functions
#' @family CAGEr clusters functions
#' 
#' @importFrom ggplot2 ggplot aes scale_fill_manual geom_histogram facet_wrap
#' @importFrom ggplot2 xlab ylab labs
#' 
#' @examples
#' plotInterquantileWidth( exampleCAGEexp, clusters = "consensusClusters"
#'                       , tpmThreshold = 50, qLow = 0.1, qUp = 0.9)
#' 
#' @export

setGeneric( "plotInterquantileWidth"
          , function( object
                    , clusters = c("tagClusters", "consensusClusters")
                    , tpmThreshold = 5
                    , qLow = 0.1, qUp = 0.9, xlim = c(0,150), ...)
            	standardGeneric("plotInterquantileWidth"))

#' @rdname plotInterquantileWidth

setMethod( "plotInterquantileWidth", "CAGEr"
         , function (object, clusters, tpmThreshold, qLow, qUp, xlim, ...) {
           
	clusters <- match.arg(clusters)
	getClustFun <- switch( clusters
	                     , tagClusters       = tagClustersGR
	                     , consensusClusters = consensusClustersGR)
  
  # Extract a list of data frames in "long" format for ggplot
	iqwidths <- lapply(seq_along(sampleLabels(object)), function(x) {
    gr <- getClustFun(object, x, returnInterquantileWidth = TRUE, qLow = qLow, qUp = qUp)
    gr <- gr[score(gr) >= tpmThreshold]
    data.frame(
      sampleName  = sampleLabels(object)[[x]],
      iq_width    = decode(gr$interquantile_width))
	})
	
	# Bind them together and set factor labels in proper order
	iqwidths <- do.call(rbind, iqwidths)
  iqwidths$sampleName <- factor(iqwidths$sampleName, levels = sampleLabels(object))
  
	binsize <- round(max(iqwidths$iq_width)/2)
	
	ggplot2::ggplot(iqwidths) +
	  ggplot2::aes(iq_width, fill = sampleName) +
	  ggplot2::scale_fill_manual(values = names(sampleLabels(object))) +
	  ggplot2::geom_histogram(bins = binsize) +
	  ggplot2::facet_wrap("~sampleName") +
	  ggplot2::xlab(paste0(
	      switch(clusters, tagClusters = "Tag Clusters", consensusClusters = "Consenss Clusters"),
	      " interquantile width (q", qLow, "-q", qUp, " (bp)")) +
	  ggplot2::ylab("Relative frequency") +
	  ggplot2::labs(fill = "Sample name")
})

#' @name plotExpressionProfiles
#' 
#' @title Plot CAGE expression profiles
#' 
#' @description Beanplot of distribution of normalized expression across CAGE
#' experiments for individual _expression classes_, colored and labeled
#' according to the information set when expression clustering was performed.
#' 
#' @param object A [`CAGEr`] object.
#' 
#' @param what `CTSS` or `consensusClusters`.
#' 
#' @details The beanplots are shown in one labeled box per _expression class_.
#' Each beanplot represents one CAGE experiment.  The vertical axis represents
#' scaled normalized expression.  The color of each class is determined by the
#' labels returned by expression clustering.
#' 
#' @author Vanja Haberle
#' @author Charles Plessy
#' 
#' @family CAGEr plot functions
#' @family CAGEr expression clustering functions
#'          
#' @examples
#' plotExpressionProfiles(exampleCAGEexp, what = "CTSS")
#' 
#' @importFrom reshape2 melt
#' @export

setGeneric( "plotExpressionProfiles"
          , function(object, what) standardGeneric("plotExpressionProfiles"))

#' @rdname plotExpressionProfiles

setMethod( "plotExpressionProfiles", "CAGEexp"
         , function (object, what=c("CTSS", "consensusClusters")) {
  what <- match.arg(what)
  
  if(what == "CTSS") {
    cl <- expressionClasses(CTSScoordinatesGR(object))
    cl.method <- CTSSexpressionClusteringMethod(object)
    DF <- CTSSnormalizedTpmDF(object)
	} else if (what == "consensusClusters") {
    cl <- expressionClasses(consensusClustersGR(object))
    cl.method <- consensusClustersExpressionClusteringMethod(object)
    DF <- DataFrame(consensusClustersTpm(object))
	}
  
  DF <- cbind(DF, exprClass=cl)
  DF <- DF[!is.na(cl),]
  df <- reshape2::melt(as.data.frame(DF), id="exprClass")

  ggplot2::ggplot(df) +
    ggplot2::aes(value, variable) +
    ggplot2::geom_violin() +
    ggplot2::facet_wrap(~exprClass) +
    ggplot2::scale_x_log10()
})


#' @name exportToBed
#' 
#' @title Create BED tracks of TSSs and clusters of TSSs
#' 
#' @description Creates BED file(s) with track(s) of individual CTSSs, tag
#' clusters or consensus clusters.  CTSSs and consensus clusters can be
#' optionally colored in the color of their expression class.  _Tag clusters_
#' and _consensus clusters_ can be displayed in a gene-like representation with
#' a line showing full span on the cluster, filled block showing interquantile
#' range and a thick box denoting position of the dominant (most frequently
#' used TSS.
#' 
#' @param object A [`CAGEr`] object.
#' 
#' @param what Which elements should be exported to BED track.  `CTSS` to export
#' individual CTSSs,  `tagClusters` to export tag clusters or `consensusClusters`
#' to export consensus clusters.
#' 
#' @param qLow,qUp Position of which "lower" (resp. "upper") quantile should be
#' used as 5' (resp. 3') boundary of the filled block in gene-like
#' representation of the cluster.  Default value `NULL` uses start (resp. end)
#' position of the cluster.  Ignored when `what = "CTSS"`.
#' 
#' @param colorByExpressionProfile Logical, should blocks be colored in the
#' color of their corresponding expression class.  Ignored when
#' `what = "tagClusters"`.
#' 
#' @param oneFile Logical, should all CAGE datasets be exported as individual
#' tracks into the same BED file (`TRUE`) or into separate BED files (`FALSE`).
#' Ignored when `what = "CTSS"`, which by default produces only one track.
#' 
#' @details The BED representations of _CTSSs_, _tag cluster_ and
#' _consensus clusters_ can be directly visualised in the ZENBU or UCSC Genome
#' Browsers.
#' 
#' When `what = "CTSS"`, one BED file with single track of 1 bp blocks
#' representing all detected CTSSs (in all CAGE samples) is created.  CTSSs can
#' be colored according to their expression class (provided the expression
#' profiling of CTSSs was done by calling [`getExpressionProfiles`] function).
#' Colors of expression classes match the colors in which they are shown in the
#' plot returned by the [`plotExpressionProfiles`] function.  For
#' `colorByExpressionProfile = FALSE`, CTSSs included in the clusters are
#' shown in black and CTSSs that were filtered out in gray.
#' 
#' When `what = "tagClusters"`, one track per CAGE dataset is created, which can
#' be exported to a single BED file (by setting `oneFile = TRUE`) or separate
#' BED files (`FALSE`).  If no quantile boundaries were provided (`qLow` and
#' `qUp` are `NULL`, TCs are represented as simple blocks showing the full
#' span of TC fromthe start to the end.  Setting `qLow` and/or `qUp` parameters
#' to a value of the desired quantile creates a gene-like representation with a
#' line showing full span of the TC, filled block showing specified
#' interquantile range and a thick 1 bp block denoting position of the dominant
#' (most frequently used) TSS.  All TCs in one track (one CAGE dataset) are
#' shown in the same color.
#' 
#' When `what = "consensusClusters"` _consensus clusters_ are exported to BED
#' file.  Since there is only one set of consensus clusters common to all CAGE
#' datasets, only one track is created in case of a simple representation.  This
#' means that when `qLow = NULL` and `qUp = NULL` one track with blocks showing
#' the full span of consensus cluster from the start to the end is created.
#' However, the distribution of the CAGE signal within consensus cluster can be
#' different in different CAGE samples, resulting in different positions of 
#' quantiles and dominant TSS.  Thus, when `qLow` and/or `qUp` parameters
#' are set to a value of the desired quantile, a separate track with a gene-like
#' representation is created for every CAGE dataset.  These tracks can be
#' exported to a single BED file (by setting `oneFile = TRUE`) or separate
#' BED files (by setting `oneFile = FALSE`).  The gene-like representation is
#' analogous to the one described above for the TCs.  In all cases consensus
#' clusters can be colored according to their expression class (provided the
#' expression profiling of consensus clusters was done by calling
#' `getExpressionProfiles` function).  Colors of expression classes match the
#' colors in which they are shown in the plot returned by the
#' `plotExpressionProfiles` function.  For `colorByExpressionProfile = FALSE`
#' all consensus clusters are shown in black.
#' 
#' @return Creates BED file(s) in the working directory.
#' 
#' @author Vanja Haberle
#' @author Charles Plessy
#' 
#' @family CAGEr export functions
#' 
#' @examples 
#' 
#' ### exporting CTSSs colored by expression class
#' exportToBed(object = exampleCAGEexp, what = "CTSS", colorByExpressionProfile = TRUE)
#' 
#' ### exporting tag clusters in gene-like representation
#'            
#' exportToBed( object = exampleCAGEexp, what = "tagClusters"
#'            , qLow = 0.1, qUp = 0.9, oneFile = TRUE)
#' 
#' @export

setGeneric( "exportToBed"
          , function( object, what = c("CTSS", "tagClusters", "consensusClusters")
                    , qLow = NULL, qUp = NULL
                    , colorByExpressionProfile = FALSE, oneFile = TRUE)
            	standardGeneric("exportToBed"))

#' @rdname exportToBed

setMethod("exportToBed", "CAGEr", function( object, what, qLow, qUp
                                          , colorByExpressionProfile, oneFile) {
	sample.labels <- sampleLabels(object)
	what <- match.arg(what)

	if(what == "CTSS") {
		
		oneFile <- TRUE
		use.blocks <- F
		ctss <- CTSScoordinatesGR(object) |> granges() |> as.data.frame()
		ctss <- data.frame(chr=ctss$seqnames, pos=ctss$start, strand=ctss$strand)
		#filtered_ctss <- object@filteredCTSSidx

		if(colorByExpressionProfile == TRUE){
			cl <- CTSSexpressionClasses(object)
			n <- names(cl)
			cl <- .extract.cluster.info(cl = cl)
			cl <- data.frame(ctss = n, x.cor = cl[,1], y.cor = cl[,2])		
			ctss <- merge(cl, ctss, by.x = "ctss", by.y = 0, all.x = T, all.y = F)
			ctss <- data.frame(ctss = ctss$ctss, chr = ctss$chr, start = ctss$pos-1, end = ctss$pos, strand = ctss$strand, x.cor = ctss$x.cor, y.cor = ctss$y.cor)
			track.file <- "CTSS.colored.by.expression.profile.bed"
			track.names <- list("CTSS (colored by expression profile)")
			clustering.method <- CTSSexpressionClusteringMethod(object)
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
		
		clusters <- tagClustersGR(object)
		
		if(!is.null(qLow)) {
		
		use.blocks <- T
		track.names <- paste(sample.labels, paste(" (tag clusters (TC) q(", qLow, ")-q(",qUp,"))", sep = ""), sep = "")
		r <- paste(".qLow", qLow, "_qUp", qUp, sep = "")
			
		}else{
			use.blocks <- F
			track.names <- paste(sample.labels, paste(" (tag clusters (TC))", sep = ""), sep = "")
			r <- ""
		}
		
		itemRgb = FALSE
		cols <- names(sample.labels)
		cols <- as.list(apply(sapply(cols, function(x) {as.integer(col2rgb(x))}), 2, function(y) {paste(y, collapse = ",")}))
		names(cols) <- sample.labels

		if(oneFile){
			track.file <- rep(paste("All.samples.tagClusters", r, ".bed", sep = ""), length(clusters))
		}else{
			track.file <- paste(sample.labels, ".tagClusters", r, ".bed", sep = "")
		}
		
	}else if(what == "consensusClusters"){
	  # consensusClustersGR not actually supported but below
	  # consensusClustersQuantileLow is also broken anyway so let's fix them together later.
		clusters <- consensusClustersGR(object)  
		colnames(clusters)[1] = "cluster"
		if(!(colorByExpressionProfile)){
			cols <- as.list(rep("0,0,0", length(sample.labels)))
		}else{
			clustering.method <- consensusClustersExpressionClusteringMethod(object)		
		cl <- consensusClustersExpressionClasses(object)
		n <- names(cl)
		cl <- .extract.cluster.info(cl = cl)
		cl <- data.frame(cluster = n, x.cor = cl[,1], y.cor = cl[,2])		
		clusters <- merge(clusters, cl, by.x = "cluster", by.y = "cluster")
		}
		
		if(!is.null(qLow) & (length(consensusClustersQuantileLow(object))>0 & length(consensusClustersQuantileUp(object))>0)) {

		if(paste("q_", qLow, sep = "") %in% colnames(consensusClustersQuantileLow(object)[[1]]) & paste("q_", qUp, sep = "") %in% colnames(consensusClustersQuantileUp(object)[[1]])){
		
		use.blocks <- T
		q.low <- consensusClustersQuantileLow(object)
		q.low <- lapply(q.low, function(x) {colnames(x)[2:ncol(x)] <- paste("qLow_", do.call(rbind, strsplit(colnames(x)[2:ncol(x)], split = "_", fixed = T))[,2], sep = ""); return(x)})
		q.up <- consensusClustersQuantileUp(object)
		q.up <- lapply(q.up, function(x) {colnames(x)[2:ncol(x)] <- paste("qUp_", do.call(rbind, strsplit(colnames(x)[2:ncol(x)], split = "_", fixed = T))[,2], sep = ""); return(x)})
		q <- lapply(as.list(1:length(sample.labels)), function(x) {merge(q.low[[x]], q.up[[x]], by.x = "cluster", by.y = "cluster")})
		names(q) <- sample.labels
						
		cumsums <- CTSScumulativesCC(object)
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
		}else if(!is.null(qLow)){
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
		names <- as.list(rep(".", length(clusters)))
	}
	
	track.descriptions <- track.names
	
	for(i in seq_along(clusters)){
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
		.make.cluster.bed.track(clusters.q = clusters[[i]], use.blocks = use.blocks, q.low = qLow, q.up = qUp, track.file = track.file[i], track.name = track.names[i], track.description = track.descriptions[i], cols = cols[[i]], name = names[[i]], itemRgb = itemRgb, app = app)
		
	}
	
	if(oneFile){
		message("\nFile '", track.file[1], "' has been created in your working directory (", getwd(), ")")	
	}else{
		message("\nFiles '", sub(sample.labels[1], "*", track.file[1]), "' for all samples have been created in your working directory (", getwd(), ")")
	}
})