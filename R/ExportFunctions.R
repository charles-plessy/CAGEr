###########################################################
# Functions for exporting results (graphical and textual)
#

.plotReverseCumulative <- function( values, col = "darkblue"
                                  , title = "", col.title = "black"
                                  , xlim = c(1, 1e5), ylim = c(1, 1e6)
                                  , xlab = "number of CAGE tags"
                                  , ylab = "number of CTSSs (>= nr tags)"
                                  , cex.axis = 1.8
                                  , add = FALSE) {
  
  # See benchmarks/sorted-abundance-benchmark.Rmd in the CAGEr's Git repository
  values <- sort(Rle(values), decreasing = TRUE)
  values <- values[values != 0]

  x <- runValue(values)
  y <- cumsum(runLength(values))
  
  if(add){
    lines(x, y, col = col, lwd = 2, type = "S")
  }else{
    plot( x, y
        , xaxt = 'n', yaxt = 'n'
        , log = "xy", type = "S"
        , lwd = 2
        , col = col
        , xlab = xlab, ylab = ylab
        , main = title
        , col.main = col.title
        , cex.axis = cex.axis, cex.lab = 1.8, cex.main = 2.5
        , xlim = xlim, ylim = ylim)	
    
    ticks <- as.integer(sapply(10^(seq(0,6,1)), function(x) {seq(x,10*x-1,x)}))
    
    axis(side = 1, at = ticks, labels = rep("", 63), cex.axis = cex.axis)   		   
    axis(side = 1, at = 10^(seq(0,6,1)), labels = formatC(10^(seq(0,6,1)), format = "f", digits = 0), cex.axis = cex.axis)		   
    axis(side = 2, at = ticks, labels = rep("", 63), cex.axis = cex.axis)   		   
    axis(side = 2, at = 10^(seq(0,6,1)), labels = formatC(10^(seq(0,6,1)), format = "f", digits = 0), cex.axis = cex.axis)
  }
}