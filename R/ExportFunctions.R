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

#' .export.bedgraph
#' @noRd
#' @importFrom rtracklayer export

.export.bedgraph <- function(data.rd, name, description, file_name, append = F) {
	data.ucsc <- as(data.rd, "UCSCData")
	data.ucsc@trackLine <- new("BasicTrackLine", name = name, description = description, visibility="full")
	score(data.ucsc) <- as.numeric(score(data.ucsc)) # Does not support Rle
	export(data.ucsc, con = file_name, format = "ucsc", subformat = "bedGraph", append = append)
}

#' @name se2grPlusOrMinus
#' 
#' @details Private function used in the export functions below.
#' 
#' @param strand "+" or "-"
#' @param data A RangedSummarisedExperiment with a single sample.
#' 
#' @return Returns a GRanges object of the required strand, where the score
#' is the expression value of the first sample of the SummarizedExperiemnt,
#' and is negative if on the negative strand.  Positions with a zero score
#' are removed.
#' 
#' @author Vanja Haberle, Charles Plessy
#' 
#' @noRd

se2grPlusOrMinus <- function(strand, data) {
  gr <- rowRanges(data)
  score(gr) <- assay(data[,1])[[1]]
  gr <- gr[score(gr) > 0]
  gr <- gr[strand(gr) == strand]
  if(strand == "-")
    score(gr) <- -1 * score(gr)
  gr
}

#' @name .export.bw.all
#' 
#' @param data A RangedSummarizedExperiment object.
#' @param sample.labels Sample labels.
#' @param v "normalized" or "raw".
#' @param genome A BSgenome object.
#' 
#' @noRd
#' @importFrom rtracklayer export.bw
#' @importFrom utils write.table

.export.bw.all <- function(data, sample.labels, v, genome) {
    rd.list <- lapply(as.list(sample.labels), function(x) {
        lapply(list("+", "-"), se2grPlusOrMinus, data[,x])
    }
    )
    names(rd.list) <- sample.labels

    a <- lapply(sample.labels, function(x){
            if(length(rd.list[[x]][[1]]) > 0){
                score(rd.list[[x]][[1]]) <- as.numeric(score(rd.list[[x]][[1]])) # Does not support Rle
                export.bw(rd.list[[x]][[1]], con = paste(x, ".CTSS.", v, ".plus.bw", sep = ""))
            }
            if(length(rd.list[[x]][[2]]) > 0){
                score(rd.list[[x]][[2]]) <- as.numeric(score(rd.list[[x]][[2]])) # Does not support Rle
                export.bw(rd.list[[x]][[2]], con = paste(x, ".CTSS.", v, ".minus.bw", sep = ""))
            }
        })
    
    description.lines = data.frame(description = unlist(lapply(sample.labels, function(x) {paste('track type=bigWig name="', paste(x, "CTSS", v, c('plus"', 'minus"'), sep = " "), ' description="', paste(x, " CTSS ", v, ' (', c("plus", "minus"), ' strand)"', sep = ""), " bigDataUrl=", paste(x, ".CTSS.", v, c(".plus", ".minus"), ".bw", sep = ""), sep = "")})))
    write.table(description.lines, file = paste("CTSS.", v, ".all.samples.track.description.txt", sep = ""), col.names = F, row.names = F, sep = "\t", quote = F)
}

.export.bedgraph.all <- function(data, sample.labels, v, oneFile) {
	rd.list <- lapply(as.list(sample.labels), function(x) {
					  lapply(as.list(c("+", "-")), se2grPlusOrMinus, data[,x])
					  }
					  )
	names(rd.list) <- sample.labels
	
	
	if(oneFile){
		if(file.exists(paste("All.samples.CTSS.", v, ".bedGraph", sep = ""))){
			r <- file.remove(paste("All.samples.CTSS.", v, ".bedGraph", sep = ""))
		}
		strands = c("plus", "minus")
		for(i in 1:length(sample.labels)){
			for(s in c(1,2)){
				if(length(rd.list[[i]][[s]])>0){
					.export.bedgraph(rd.list[[i]][[s]], name = paste(sample.labels[i], "_", v, "_", strands[s], sep = ""), description = paste(sample.labels[i], " CTSS ", v, " (", strands[s], " strand)", sep = ""), file_name = paste("All.samples.CTSS.", v, ".bedGraph", sep = ""), append = T)
				}
			}
		}
	}else{
		a <- lapply(as.list(sample.labels), function(x) {
				if(length(rd.list[[x]][[1]]) > 0){
					.export.bedgraph(rd.list[[x]][[1]], name = paste(x, "_", v, "_plus", sep = ""), description = paste(x, " CTSS ", v, " (plus strand)", sep = ""), file_name = paste(x, ".CTSS.", v, ".plus.bedGraph", sep = ""), append = F)
				}
				if(length(rd.list[[x]][[2]]) > 0){	
					.export.bedgraph(rd.list[[x]][[2]], name = paste(x, "_", v, "_minus", sep = ""), description = paste(x, " CTSS ", v, " (minus strand)", sep = ""), file_name = paste(x, ".CTSS.", v, ".minus.bedGraph", sep = ""), append = F)			   
				}
			   }
			   )
	}
		
}
