#' @include CAGEr.R

#' @name hanabi-class
#' 
#' @description TBD
#' 
#' @title Hanabi class
#' 
#' @details TBD
#' 
#' @rdname hanabi-class
#' 
#' @import methods
#' @importFrom S4Vectors List

.hanabi <- setClass("hanabi", contains = "SimpleList", prototype = c())
# (See <https://github.com/Bioconductor/Contributions/issues/261#issuecomment-277479436>.)

#' @name hanabi
#' 
#' @title Calcultate richness in preparation for plotting
#' 
#' @description Rarefy data at multiple sample sizes using the
#' \code{vegan} package and return a \sQuote{hanabi} object that can be passed
#' to plot functions.
#' 
#' The computation can be long, so the steps of rarefaction and plotting are kept
#' separate.
#' 
#' @param x An object contained expression counts on which richness scores can
#'        be calculated.  For example an expression table in \code{DataFrame}
#'        or \code{data.frame} format where columns are samples and rows are
#'        featuressuch as genes, TSS, etc, or a vector of counts (tag counts,
#'        molecule counts, ...), or \code{GRanges} or \code{GRangesList}
#'        objects, etc.
#'        
#' @param n The maximum number of rarefactions per sample.
#' 
#' @param step Subsample sizes are calculated by taking the largest sample
#'        and multiplying it by the step "n" times.
#'        
#' @param from Add one sample size (typically "0") in order to extend the
#'        plot on the left-hand side.
#'        
#' @param useMulticore Logical, should multicore be used.
#'   \code{useMulticore = TRUE} has no effect on non-Unix-like platforms.  At
#'   the moment, it also has only effects on lists and list-derived classes
#'   (data frames but not matrices).
#'   
#' @param nrCores Number of cores to use when \code{useMulticore = TRUE}
#'   (set to \code{NULL} to use all detected cores).
#' 
#' @details This function does not take directly CAGEr objects as input,
#' because hanabi plots can be made from CTSS, clustered or gene-level
#' data, therefore it is not possible to guess which one to use.
#'
#' @return A list-based object of class "hanabi".
#'
#' @family CAGEr richness functions
#' @seealso \code{vegan::rarecurve}.
#' 
#' @author Charles Plessy
#' 
#' @importFrom vegan rarefy
#' @importFrom S4Vectors DataFrame Rle
#' @export hanabi
#' 
#' @examples
#' h <- hanabi(CTSStagCountDF(exampleCAGEexp))
#' h
#' plot(h)
#' hanabi(CTSStagCountGR(exampleCAGEexp, 2))

setGeneric( "hanabi"
          , function( x
                    , n    = 20
                    , step = 0.75
                    , from = NULL
                    , useMulticore = FALSE
                    , nrCores      = NULL)
              standardGeneric("hanabi"))

#' @rdname hanabi

setMethod("hanabi", "Rle", function(x, n, step, from, useMulticore, nrCores) {
  hanabi(as.vector(x), n = n, step = step, from = from)
})

#' @rdname hanabi

setMethod("hanabi", "numeric", function(x, n, step, from, useMulticore, nrCores) {
  if (any(round(x) != x))
    stop("All scores must be integers.")
  hanabi(as.integer(x))
})

#' @rdname hanabi

setMethod("hanabi", "integer", function(x, n, step, from, useMulticore, nrCores) {
  ns <- step ^ (0:n)
  ns <- round(sum(x) * ns)
  if (! is.null(from))
    ns <- c( ns[ns > from], from)
  ntags <- sum(x)
  ns <- c(ntags, ns[ns < ntags])
  h <- SimpleList()
  h[[1]] <- xy.coords(x = ns, y = rarefy(x, ns))
  .hanabi(h)
})

#' @rdname hanabi

setMethod("hanabi", "GRanges", function(x, n, step, from, useMulticore, nrCores) {
  if (any(round(score(x)) != score(x)))
    stop("All scores must be integers.")
    hanabi(score(x))
})

# Same method for list-like or List-like objects, such as list,
# data.frame, DataFrame or GRangesList.
# Uses `unlist(SimpleList())` to move the hanabi class from inside to
# outside the list, taking advantage of the base class SimpleList.

#' @rdname hanabi

setMethod("hanabi", "List", function(x, n, step, from, useMulticore, nrCores) {
  l <- bplapply( x, hanabi, n = n, step = step, from = from
               , BPPARAM = CAGEr_Multicore(useMulticore, nrCores))
  unlist(SimpleList(l))
})

#' @rdname hanabi

setMethod("hanabi", "list", function(x, n, step, from, useMulticore, nrCores) {
  l <- bplapply( x, hanabi, n = n, step = step, from = from
               , BPPARAM = CAGEr_Multicore(useMulticore, nrCores))
  unlist(SimpleList(l))
})

#' @rdname hanabi

setMethod("hanabi", "matrix", function(x, n, step, from, useMulticore, nrCores) {
  unlist(SimpleList(apply(x, 2, hanabi, n = n, step = step, from = from)))
})


#' @name .add.alpha
#' @noRd
#' @description Accessory function to make the lines a little transparent.
#' @seealso https://gist.github.com/mages/5339689#file-add-alpha-r
#' @importFrom grDevices col2rgb

.add.alpha <- function(col, alpha) {
  apply( sapply(col, col2rgb) / 255
       , 2
       , function(x)
           rgb(x[1], x[2], x[3], alpha=alpha))
}

#' @name plot.hanabi
#' 
#' @title Plotting Hanabi objects
#' 
#' @description S3 method to plot hanabi objects.  Used by the
#' \code{\link{hanabiPlot}} function.
#' 
#' @param x The hanabi object to plot.
#' @param alpha The alpha transparency of the plot lines.
#' @param col A vector indicating a color per sample (or a vector that
#'        can be recycled that way).
#' @param xlab Horizontal axis label.
#' @param ylab Vertical axis label.
#' @param main Plot title.
#' @param pch Plot character at the tip of the lines.
#' @param ... Other parameters passed to the generic plot, points or lines functions.
#' 
#' @family CAGEr richness functions
#' 
#' @author Charles Plessy
#' 
#' @importFrom graphics lines
#' 
#' @export

plot.hanabi <-
  function( x
          , alpha = 0.5
          , col   = "black"
          , xlab  = "Total counts"
          , ylab  = "Unique features"
          , main  = "Hanabi plot"
          , pch   = 1
          , ...) {
  xmax <- sapply(x, function(x) max(x$x))
  xmin <- sapply(x, function(x) min(x$x))
  ymax <- sapply(x, function(x) max(x$y))
  ymin <- sapply(x, function(x) min(x$y))
  plot( c(min(xmin), max(xmax))
      , c(min(ymin), max(ymax))
      , type="n"
      , xlab = xlab
      , ylab = ylab
      , main = main
      , ...)
  lines( x
       , col = .add.alpha(col, alpha))
  points( x
        , col = col
        , pch = pch)
}

#' @name points.hanabi
#' @rdname plot.hanabi
#' @export

points.hanabi <- function(x, ...) {
  xmax <- sapply(x, function(x) max(x$x))
  ymax <- sapply(x, function(x) max(x$y))
  points(xmax, ymax, ...)
}

#' @name lines.hanabi
#' @rdname plot.hanabi
#' @export

lines.hanabi  <- function(x, ...) {
  invisible(Map(lines, x, ...))
}

#' hanabiPlot
#' 
#' Plot feature discovery curves
#' 
#' Plots the number of features (genes, transcripts, ...) detected for a
#' given number of counts (reads, unique molecules, ...).  Each library is
#' sub-sampled by rarefaction at various sample sizes, picked to provide
#' enough points so that the curves look smooth.  The final point is plotted
#' as an open circle, hence the name "hanabi", which means fireworks in
#' Japanese. 
#' 
#' The rarefactions take time to do, so this step is done by a separate
#' function, so that the result is easily cached.
#' 
#' @param x A hanabi object.
#' @param group A character vector or a factor grouping the samples.
#' @param col A character vector colors (at most one per group).
#' @param legend.pos Position of the legend, passed to the \code{legend} function.
#' @param pch Plot character at the tip of the lines and in the legend.
#' @param ... Further arguments to be passed to the \code{plot.hanabi} function.
#' 
#' @family CAGEr richness functions
#' @author Charles Plessy
#'  
#' @examples
#' h <- hanabi(CTSStagCountDF(exampleCAGEexp))
#' hanabiPlot(h, group = 1:5)
#' hanabiPlot(hanabi(CTSStagCountDF(exampleCAGEexp), n = 20, step = 0.8, from = 25000), group = 1:5)
#' hanabiPlot(hanabi(CTSStagCountDF(exampleCAGEexp), n = 10, step = 0.98), group = 1:5)
#' hanabiPlot(h, group=c("A", "A", "B", "C", "B"), col=c("red", "green", "blue"))
#' hanabiPlot(h, group = 1:5, pch=1:5, col="purple")
#' 
#' @family CAGEr richness functions
#' @family CAGEr plot functions
#' 
#' @export hanabiPlot

hanabiPlot <- function ( x
                       , group
                       , col = NULL
                       , legend.pos = "topleft"
                       , pch = 1
                       , ...) {
  if(missing(group))
    stop("Missing ", dQuote("group"), " argument. If no group is needed use directly the ",
         dQuote("plot()"), " function instead.")
  
  # When coercing to factor, make sure that 1) levels are in their order of
  # appearance in the original character vector and 2) levels are not reordered
  # if the coerced object was already a factor.
  if (! is.factor(group))
    group <- factor(group, unique(group))
  
  if (length(col) > nlevels(group))
    stop("More colors than levels!")
  
  if (is.null(col)) {
    # Take group levels as colors if no colors were provided.
    col <- group
    levels(col) <- 1:nlevels(col)
  } else {
    # Recycle color vector and give color to the groups.
    col.r <- rep(col, length.out = nlevels(group))
    col <- group
    levels(col) <- col.r
  }

  plot(x, pch = pch, col = as.character(col), ...)
  legend( x = legend.pos
          , legend = levels(group)
          , col = levels(col)
          , pch = pch)
  invisible()
}
