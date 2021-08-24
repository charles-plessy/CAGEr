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
#' @param xlim Range of width to be plotted.
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
#' @importFrom ggplot2 ggplot scale_fill_manual geom_histogram facet_wrap
#' @importFrom ggplot2 ggtitle xlab ylab labs
#' 
#' @examples
#' 
#' plotInterquantileWidth( exampleCAGEexp, clusters = "tagClusters"
#'                       , tpmThreshold = 50, qLow = 0.1, qUp = 0.9
#'                       , xlim = c(2,200))
#'                       
#' plotInterquantileWidth( exampleCAGEexp, clusters = "consensusClusters"
#'                       , tpmThreshold = 50, qLow = 0.1, qUp = 0.9
#'                       , xlim = c(2,200))
#' 
#' @export

setGeneric( "plotInterquantileWidth"
          , function( object
                    , clusters = c("tagClusters", "consensusClusters")
                    , tpmThreshold = 5
                    , qLow = 0.1, qUp = 0.9, xlim = c(0,150))
            	standardGeneric("plotInterquantileWidth"))

#' @rdname plotInterquantileWidth

setMethod( "plotInterquantileWidth", "CAGEexp"
         , function (object, clusters, tpmThreshold, qLow, qUp, xlim) {
           
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
  iqwidths <- iqwidths[iqwidths$iq_width >= xlim[1] & iqwidths$iq_width <= xlim[2],]
  
	binsize <- round(max(iqwidths$iq_width)/2)
	
	ggplot2::ggplot(iqwidths) +
	  ggplot2::aes_string(x = "iq_width", fill = "sampleName") +
	  ggplot2::scale_fill_manual(values = names(sampleLabels(object))) +
	  ggplot2::geom_histogram(bins = binsize) +
	  ggplot2::facet_wrap("~sampleName") +
	  ggplot2::ggtitle(paste0(
	    switch(clusters, tagClusters = "Tag Clusters", consensusClusters = "Consenss Clusters"),
	    " interquantile width (quantile ", qLow, " to ", qUp, ")")) +
	  ggplot2::xlab("Interquantile width (bp)") +
	  ggplot2::ylab("Frequency") +
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
#' exampleCAGEexp |> plotExpressionProfiles("consensusClusters")
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
    ggplot2::aes_string(x = "value", y = "variable") +
    ggplot2::geom_violin() +
    ggplot2::facet_wrap(~exprClass) +
    ggplot2::scale_x_log10()
})


#' @name exportToTrack
#' 
#' @title Converts TSSs and clusters of TSSs to a genome browser track format
#' 
#' @description Converts _CTSS_, _tag clusters_ or _consensus clusters_ to the
#' `UCSCData` format of the `rtracklayer` package, that can be exported to BED
#' file(s) with track information for genome browsers.  _CTSSes_ and
#' _consensus clusters_ are optionally colored by their expression class.
#' _Tag clusters_ and _consensus clusters_ can be displayed in a
#' whiskerplot-like representation with a line showing full span on the cluster,
#' filled block showing interquantile range and a thick box denoting position of
#' the dominant (most frequently) used _TSS_.
#' 
#' @param object A [`CAGEexp`] object.
#' 
#' @param what Which elements should be exported: `CTSS` for
#' individual _CTSSs_,  `tagClusters` for _tag clusters_ or `consensusClusters`
#' for _consensus clusters_.
#' 
#' @param qLow,qUp Position of which "lower" (resp. "upper") quantile should be
#' used as 5' (resp. 3') boundary of the filled block in whiskerplot-like
#' representation of the cluster.  Default: `NULL` (plain line representation).
#' Ignored when `what = "CTSS"`.
#' 
#' @param colorByExpressionProfile Logical, should blocks be colored in the
#' color of their corresponding expression class.  Ignored when
#' `what` equals`"tagClusters"`.
#' 
#' @param oneTrack Logical, should the data be converted in an individual
#' object or a list of objects?
#' 
#' @details The BED representations of _CTSSs_, _tag cluster_ and
#' _consensus clusters_ can be directly visualised in the ZENBU or UCSC Genome
#' Browsers.
#' 
#' When `what = "CTSS"`, one `UCSCData` object with single track of 1 bp blocks
#' representing all detected CTSSs (in all CAGE samples) is created.  CTSSs can
#' be colored according to their expression class (see
#' [`getExpressionProfiles`] and [`plotExpressionProfiles`]).  For
#' `colorByExpressionProfile = FALSE`, CTSSs included in the clusters are
#' shown in black and CTSSs that were filtered out in gray.
#' 
#' When `what = "tagClusters"`, one track per CAGE dataset is created, which can
#' be exported to a single `UCSCData` object (by setting `oneFile = TRUE`) or separate
#' ones (`FALSE`).  If no quantile boundaries were provided (`qLow` and
#' `qUp` are `NULL`, TCs are represented as simple blocks showing the full
#' span of TC fromthe start to the end.  Setting `qLow` and/or `qUp` parameters
#' to a value of the desired quantile creates a gene-like representation with a
#' line showing full span of the TC, filled block showing specified
#' interquantile range and a thick 1 bp block denoting position of the dominant
#' (most frequently used) TSS.  All TCs in one track (one CAGE dataset) are
#' shown in the same color.
#' 
#' When `what = "consensusClusters"` _consensus clusters_ are exported.
#' Since there is only one set of consensus clusters common to all CAGE
#' datasets, only one track is created in case of a simple representation.  This
#' means that when `qLow = NULL` and `qUp = NULL` one track with blocks showing
#' the full span of consensus cluster from the start to the end is created.
#' However, the distribution of the CAGE signal within consensus cluster can be
#' different in different CAGE samples, resulting in different positions of 
#' quantiles and dominant TSS.  Thus, when `qLow` and/or `qUp` parameters
#' are set to a value of the desired quantile, a separate track with a gene-like
#' representation is created for every CAGE dataset.  These tracks can be
#' exported to a single `UCSCData` object (by setting `oneFile = TRUE`) or separate
#' ones (by setting `oneFile = FALSE`).  The gene-like representation is
#' analogous to the one described above for the TCs.  In all cases consensus
#' clusters can be colored according to their expression class (provided the
#' expression profiling of consensus clusters was done by calling
#' `getExpressionProfiles` function).  Colors of expression classes match the
#' colors in which they are shown in the plot returned by the
#' `plotExpressionProfiles` function.  For `colorByExpressionProfile = FALSE`
#' all consensus clusters are shown in black.
#' 
#' @return Returns either a `rtracklayer` `UCSCData` object, or a `GRangesList`
#' of them.
#' 
#' @author Vanja Haberle
#' @author Charles Plessy
#' 
#' @family CAGEr export functions
#' 
#' @examples 
#' # You can export from a CAGEexp object or from a cluster object directly:
#' exportToTrack(exampleCAGEexp, what = "CTSS")  # Is same as:
#' exportToTrack(CTSScoordinatesGR(exampleCAGEexp))  # Or:
#' exampleCAGEexp |> CTSScoordinatesGR() |> exportToTrack()
#' 
#' # Export a single sample, 
#' exampleCAGEexp |> CTSStagCountGR(2)      |> exportToTrack()
#' exampleCAGEexp |> CTSSnormalizedTpmGR(2) |> exportToTrack()
#' 
#' # Exporting multiple samples results in a GRangesList of UCSCData objects.
#' exportToTrack(exampleCAGEexp, what = "CTSS", oneTrack = FALSE)
#' exampleCAGEexp |> CTSStagCountGR("all")  |> exportToTrack()
#' exampleCAGEexp |> CTSSnormalizedTpmGR("all")  |> exportToTrack()
#' 
#' ### exporting CTSSs colored by expression class
#' # Temporarly disabled
#' # exportToTrack(exampleCAGEexp, what = "CTSS", colorByExpressionProfile = TRUE)
#' 
#' ### exporting tag clusters in gene-like representation
#' exportToTrack(exampleCAGEexp, what = "tagClusters", qLow = 0.1, qUp = 0.9)
#' tagClustersGR(exampleCAGEexp, 1) |> exportToTrack(qLow = 0.1, qUp = 0.9)
#'            
#' ### exporting consensus clusters
#' exportToTrack( exampleCAGEexp, what = "consensusClusters")
#' exampleCAGEexp |>
#'   consensusClustersGR("Zf.high", qLow = .1, qUp = .9) |>
#'   exportToTrack(qLow = .1, qUp = .9)
#' 
#' @export

setGeneric( "exportToTrack"
          , function( object, what = c("CTSS", "tagClusters", "consensusClusters")
                    , qLow = NULL, qUp = NULL
                    , colorByExpressionProfile = FALSE, oneTrack = TRUE)
            standardGeneric("exportToTrack"))

#' @rdname exportToTrack

setMethod( "exportToTrack", "CAGEexp"
         , function( object, what, qLow, qUp
                   , colorByExpressionProfile, oneTrack) {
  clusters <-
	  switch( match.arg(what)
	        , CTSS              = if (oneTrack) { CTSScoordinatesGR(object)
	                              } else { CTSStagCountGR(object, samples = "all") }
	        , tagClusters       = tagClustersGR(      object, qLow = qLow, qUp = qUp)
	        , consensusClusters = consensusClustersGR(object, qLow = qLow, qUp = qUp))
	
	exportToTrack( clusters, qLow = qLow, qUp = qUp
	             , colorByExpressionProfile = colorByExpressionProfile
	             , oneTrack = oneTrack)
})

#' @rdname exportToTrack

setMethod( "exportToTrack", "GRangesList"
           , function( object, what, qLow, qUp
                       , colorByExpressionProfile, oneTrack) {
  grl <- endoapply( object, exportToTrack
                  , qLow = qLow, qUp = qUp
                  , colorByExpressionProfile = colorByExpressionProfile
                  , oneTrack = FALSE)
  if (isTRUE(oneTrack)) {
    unlist(grl)
  } else {
    GRangesList(grl, compress = FALSE) # Because rtracklayer::export.bed does not like compressed GRangesList...
  }
})

#' @importFrom IRanges IRangesList
#' @rdname exportToTrack

setMethod( "exportToTrack", "GRanges",
function(object, what, qLow, qUp, colorByExpressionProfile, oneTrack) {
  if((! is.null(qLow)) & (! is.null(qUp)) ) {
    decoded.mat <-
      rbind( width(object)
           , mcols(object)[,paste0("q_", qLow)] |> decode()
           , mcols(object)[,paste0("q_", qUp)]  |> decode()
      )
    object$blocks <- IRangesList(apply(decoded.mat, 2, \(x) {
      # See benchmarks/BED12-blockInfo-benchmark.md
      width_value <- x[1]
      qLow_value  <- x[2]
      qUp_value   <- x[3]
      ir <- IRanges()                      |>
        c( if(qLow_value != 1) IRanges(1)) |>
        c( IRanges(qLow_value, qUp_value)) |>
        c( if(qUp_value != width_value) IRanges(width_value))
    }))
  }
  ucsc <- as(object, "UCSCData")
  score(ucsc) <- decode(score(ucsc))
  ucsc@trackLine <- new( "BasicTrackLine", name = "TC"
                       , description = "CAGE Tag Clusters (TC)"
                       , visibility="full")
  ucsc
})

#' @rdname exportToTrack

setMethod( "exportToTrack", "CTSS",
function( object, what, qLow, qUp, colorByExpressionProfile, oneTrack) {
  if (isTRUE(colorByExpressionProfile)) {
    stop("Coloring by expression profile is disabled for the moment.")
  } else {
    object$itemRgb <- ifelse(object$filteredCTSSidx, "black", "grey50")
  }
  if (is.null(score(object))) score(object) <- 0L
  trk <- exportToTrack( GRanges(object), qLow = qLow, qUp = qUp
                      , colorByExpressionProfile = colorByExpressionProfile
                      , oneTrack = oneTrack)
  sampleLabel <- sampleLabels(object)
  if(is.null(sampleLabel)) sampleLabel <- "All samples pooled"
  trk@trackLine@name <- paste0(sampleLabel, " (", trk@trackLine@name, ")")
  trk@trackLine@description <- paste0(sampleLabel, " (", trk@trackLine@description, ")")
  trk
})

#' @rdname exportToTrack

setMethod( "exportToTrack", "TagClusters",
function( object, what, qLow, qUp, colorByExpressionProfile, oneTrack) {
  object$thick <- IRanges(object$dominant_ctss)
  object$dominant_ctss <- NULL
 	names(object) <- NULL
 	object$name <- NA
 	object$nr_ctss <- NULL
 	object$tpm.dominant_ctss <- NULL
  score(object) <- 0L
  exportToTrack( GRanges(object), qLow = qLow, qUp = qUp
               , colorByExpressionProfile = colorByExpressionProfile
               , oneTrack = oneTrack)
})

#' @rdname exportToTrack

setMethod( "exportToTrack", "ConsensusClusters"
           , function( object, what, qLow, qUp
                       , colorByExpressionProfile, oneTrack) {
  score(object) <- 0L
  exportToTrack( GRanges(object), qLow = qLow, qUp = qUp
                 , colorByExpressionProfile = colorByExpressionProfile
                 , oneTrack = oneTrack)
})

# This function might be used later to restore the capacity of assigning a
# color to each SOM cluster.

.myColorMatrix <- function(color.mx, nrows, ncols, ...) {
  .myColorRamp <- function (vec, color.low="green", color.high="red", color.mid=NULL, alpha=1, value.low=min(vec), value.high=max(vec), value.mid=(value.low+value.high)/2, ...) {
    vec.01 <- rep(NA, length(vec))
    vec.01[vec <= value.low] <- 0
    vec.01[vec >= value.high] <- 1
    vec.01[vec>value.low & vec<=value.mid] <- 0.5*(vec[vec>value.low & vec<=value.mid]-value.low)/(value.mid-value.low)
    vec.01[vec>value.mid & vec<value.high] <- 0.5+ 0.5*(vec[vec>value.mid & vec<value.high]-value.mid)/(value.high-value.mid)
      cr <- colorRamp(c(color.low, color.mid, color.high),...)
      return(apply (cr(vec.01)/255, 1, function(x) rgb(x[1], x[2], x[3], alpha)))
    }
  top.row <- .myColorRamp(1:ncols,
    color.low=color.mx[1,1],
    color.mid=NULL,
    color.high=color.mx[1,2], ...)
  bottom.row <- .myColorRamp(1:ncols,
    color.low=color.mx[2,1],
    color.mid=NULL,
    color.high=color.mx[2,2], ...)
  sapply(1:ncols, function(i) .myColorRamp(1:nrows,
    color.low=top.row[i],
    color.mid=NULL,
    color.high=bottom.row[i], ...))
}