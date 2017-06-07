#' @include AllClasses.R CAGEexp.R

#' @name hanabi-class
#' 
#' @title Hanabi class
#' 
#' @details TBD
#' 
#' @rdname CAGEexp-class
#' 
#' @import methods
#' @importFrom S4Vectors List

.hanabi <- setClass("hanabi", contains = "SimpleList", prototype = c())
# (See <https://github.com/Bioconductor/Contributions/issues/261>.)

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
#' @param x An expression table where columns are samples and rows are features
#'        such as genes, TSS, etc, or a vector of counts (tag counts, molecule
#'        counts, ...).
#'        
#' @param n The maximum number of rarefactions per sample.
#' 
#' @param step Subsample sizes are calculated by taking the largest sample
#'        and multiplying it by the step "n" times.
#'        
#' @param from Add one sample size (typically "0") in order to extend the
#'        plot on the left-hand side.
#' 
#' @details This function does not take directly CAGEr objects as input,
#' because hanabi plots can be made from CTSS, clustered or gene-level
#' data, therefore is not possible to guess which one to use.
#'
#' @return A list-based object of class "hanabi".
#'
#' @family CAGEr richness functions
#' @seealso \code{\link{vegan::rarecurve}}.
#' 
#' @author Charles Plessy
#' 
#' @importFrom vegan rarefy
#' @importFrom S4Vectors DataFrame Rle
#' @export hanabi
#' 
#' @examples
#' ce <- readRDS(system.file(package = "CAGEr", "extdata/CAGEexp.rds"))
#' h <- hanabi(CTSStagCountDF(ce))
#' h
#' plot(h)

setGeneric( "hanabi"
          , function( x
                    , n    = 20
                    , step = 0.75
                    , from = NULL)
              standardGeneric("hanabi"))

setMethod("hanabi", "Rle", function(x, n, step, from) {
  hanabi(as.vector(x), n = n, step = step, from = from)
})

setMethod("hanabi", "integer", function(x, n, step, from) {
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

# Same method for DataFrame and data.frame.
# Uses `unlist(SimpleList())` to move the hanabi class from inside to
# outside the list, taking advantage of the base class SimpleList.

setMethod("hanabi", "DataFrame", function(x, n, step, from) {
  unlist(SimpleList(lapply(d, hanabi, n = n, step = step, from = from)))
})

setMethod("hanabi", "data.frame", function(x, n, step, from) {
  unlist(SimpleList(lapply(d, hanabi, n = n, step = step, from = from)))
})

#' points.hanabi
#' 
#' Add a final point in hanabi plots.
#' 
#' Will only add a point for the final, non-subsampled value of each
#' sample of in a hanabi object.
#' 
#' @param h The hanabi object.
#' @param ... Other parameters passed to the generic points function
#' 
#' @family CAGEr richness functions
#' @author Charles Plessy
#' 
#' @export points.hanabi
#' @export lines.hanabi

points.hanabi <- function(h, ...) {
  xmax <- sapply(h, function(x) max(x$x))
  ymax <- sapply(h, function(x) max(x$y))
  points(xmax, ymax, ...)
}

lines.hanabi  <- function(h, ...) {
  invisible(Map(lines, h, ...))
}

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

#' plot.hanabi
#' 
#' Plotting Hanabi objects
#' 
#' @param h The hanabi object to plot.
#' @param alpha The alpha transparency of the plot lines.
#' @param col A vector indicating a color per sample (or a vector that
#'        can be recycled that way).
#' @param xlab Horizontal axis label.
#' @param ylab Vertical axis label.
#' @param main Plot title.
#' @param pch Plot character at the tip of the lines.
#' @param ... other arguments passed to the generic plot function.
#' 
#' @family CAGEr richness functions
#' 
#' @author Charles Plessy
#' 
#' @export plot.hanabi

plot.hanabi <-
  function( h
          , alpha = 0.5
          , col   = "black"
          , xlab  = "Total counts"
          , ylab  = "Unique features"
          , main  = "Hanabi plot"
          , pch   = 1
          , ...) {
  xmax <- sapply(h, function(x) max(x$x))
  xmin <- sapply(h, function(x) min(x$x))
  ymax <- sapply(h, function(x) max(x$y))
  ymin <- sapply(h, function(x) min(x$y))
  plot( c(min(xmin), max(xmax))
      , c(min(ymin), max(ymax))
      , type="n"
      , xlab = xlab
      , ylab = ylab
      , main = main
      , ...)
  lines( h
       , col = .add.alpha(col, alpha))
  points( h
        , col = col
        , pch = pch)
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
#' @param group A vector grouping the samples.  Coerced to factor.
#' @param col A vector of colors
#' @param legend.pos Position of the legend, passed to the \code{legend} function.
#' @param pch Plot character at the tip of the lines and in the legend.
#' @param ... Further arguments to be passed to the \code{plot.hanabi} function.
#' 
#' @family CAGEr richness functions
#' @author Charles Plessy
#'  
#' @examples
#' ce <- readRDS(system.file(package = "CAGEr", "extdata/CAGEexp.rds"))
#' h <- hanabi(CTSStagCountDF(ce))
#' hanabiPlot(h, group = 1:5)
#' hanabiPlot(hanabi(CTSStagCountDF(ce), n = 20, step = 0.8, from = 0))
#' hanabiPlot(hanabi(CTSStagCountDF(ce), n = 10, step = 0.95))
#' hanabiPlot(h, group = c("A", "A", "B", "C", "B"), col=c("red", "green", "blue"))
#' hanabiPlot(h, col="purple")
#' 
#' @family CAGEr richness functions
#' @family CAGEr plot functions
#' 
#' @importFrom vegan rarefy
#' @export hanabiPlot

hanabiPlot <- function ( x
                       , group=NULL
                       , col = "black"
                       , legend.pos = "topleft"
                       , pch
                       , ...) {
  # Coerce group to factor
  if (! is.null(group)) group <- factor(group)
  
  # Take user-provided color, or take group levels as colors.
  
  if (missing(col) & ! is.null(group)) {
    col  <- 1:nlevels(group)
    cols <- as.numeric(group)
  }

  plot(x, pch = pch, col = col, ...)
  
  if (! is.null(group)) {
    legend( x = legend.pos
          , legend = levels(group)
          , col = levels(factor(col))
          , pch = pch)
  }
  invisible()
}
