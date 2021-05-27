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



.myColorRamp <- function (vec, color.low="green", color.high="red", color.mid=NULL, alpha=1,
value.low=min(vec), value.high=max(vec), value.mid=(value.low+value.high)/2, ...) {
    vec.01 <- rep(NA, length(vec))
    vec.01[vec <= value.low] <- 0
    vec.01[vec >= value.high] <- 1
    vec.01[vec>value.low & vec<=value.mid] <-
	0.5*(vec[vec>value.low & vec<=value.mid]-value.low)/(value.mid-value.low)
    vec.01[vec>value.mid & vec<value.high] <-
	0.5+ 0.5*(vec[vec>value.mid & vec<value.high]-value.mid)/(value.high-value.mid)
    cr <- colorRamp(c(color.low, color.mid, color.high),...)
    return(apply (cr(vec.01)/255, 1, function(x) rgb(x[1], x[2], x[3], alpha)))
}

.myColorMatrix <- function(color.mx, nrows, ncols, ...) {
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
											 color.high=bottom.row[i], ...
											 
											 )
		   )
	
}



.extract.cluster.info <- function(tpm.mx = NULL, cl) {

	n <- names(cl)
	if(length(grep("_", cl, fixed = T)) > 0) {
		cl <- strsplit(cl, split = "_", fixed = T)
		cl <- do.call(rbind, cl)
		cl <- apply(cl, 2, as.integer)		
	}else{
		cl <- as.matrix(cbind(as.integer(cl)-1, rep(0, length(cl))))
	}
	
	if(!(is.null(tpm.mx))){
	m <- t(scale(t(log(tpm.mx+1)), center=F)) 
	m <- m[rownames(tpm.mx) %in% as.integer(n), ]
	
	return(list(cl, m))
	}else{
		return(cl)
	}

}


#####
# Function that plots beanplots of distribution of signal for different SOM classes of CTSSs across given samples/stages (SOM classes map) 
# ARGUMENTS: value matrix - matrix containing expression values for each CTSS (rows) across different stages/samples (columns)
#            y.som - SOM object returned by calling som clustering function on above given value matrix
#			 dim.som.x - number of classes in the x-dimension
#			 dim.som.y - number of classes in the y-dimention
#			 plot.file - path of the plot file
#			 cols - vector of 4 colors (corners: bottomleft, topleft, bottomright, topright) for designing a color matrix for SOM map
# RETURNS: plots beanplots of distribution of signal for different SOM classes of CTSSs across given samples/stages (SOM classes map)

#' @name .plot.clusters.beanplots
#' @noRd
#' @importFrom beanplot beanplot

.plot.clusters.beanplots <- function(value.matrix, cl, cl.method, dim.som.x, dim.som.y, ylim = c(0,2), las = 0, labels = colnames(value.matrix), titles = 'number', cex.axis = 1, cex.main = 2, cex.lab = 1, cols = c("red", "gold", "green", "blue")) {
	
	color.matrix.solid <- .myColorMatrix(matrix(cols, nrow=2), nrows=dim.som.y, ncols=dim.som.x)
	color.matrix.solid = matrix(as.vector(color.matrix.solid), ncol = dim.som.x, byrow = F)
	
	for (j in (dim.som.y:1)-1) {
		for (i in (1:dim.som.x)-1) {
			to.plot <- data.frame(value.matrix[cl[,1] == i & cl[,2] ==j,])
			if(titles == 'class'){
				if(dim.som.y == 1){
					title = i + 1
				}else{
					title = paste(i, j, sep=",")
				}
			}
			if(titles == 'number'){
				if(cl.method == "som"){
					title = paste(i, "_", j, " (", nrow(to.plot), ")", sep = "")
				}else if(cl.method == "kmeans"){
					title = paste(i, " (", nrow(to.plot), ")", sep = "")
				}
			}
			if (nrow(to.plot)>2000) to.plot <- to.plot[sample(nrow(to.plot), 2000),]
			if(j==0){
				beanplot(to.plot, ylim=ylim, bw=0.1, las=las, log = "", beanlinewd=0, cex.axis=cex.axis, col=c(color.matrix.solid[j+1,i+1],rgb(0,0,0,0),rgb(0,0,0,0),'white'), border=NA, main=title, yaxt = 'n', names = labels, cex.lab = cex.lab, cex.main = cex.main, col.main = color.matrix.solid[j+1,i+1])
			}else{
				beanplot(to.plot, ylim=ylim, bw=0.1, las=las, log = "", beanlinewd=0, cex.axis=cex.axis, col=c(color.matrix.solid[j+1,i+1],rgb(0,0,0,0),rgb(0,0,0,0),'white'), border=NA, main=title, yaxt = 'n', show.names = F, cex.lab = cex.lab, cex.main = cex.main, col.main = color.matrix.solid[j+1,i+1])
				
			}
			
		}
	}
	
}

#' @importFrom IRanges reverse

.make.cluster.bed.track <- function(clusters.q, use.blocks = T, q.low = NULL, q.up = NULL, track.file, track.name = 'q_track', track.description = 'q track', name = ".", cols = "0,0,0", itemRgb = FALSE, app = T) {
	
	if(!(app) & file.exists(track.file)){
		file.remove(track.file)
	}	

	chr = seqnames(clusters.q)
	start = start(clusters.q)
	end = end(clusters.q)
	strand = strand(clusters.q)

	if(use.blocks){
		
	q.low.pos = decode(mcols(clusters.q)[,paste('q_', q.low, sep = '')] - 1 )
	q.up.pos  = decode(mcols(clusters.q)[,paste('q_', q.up,  sep = '')]     )
		
	dominant_ctss_start = mcols(clusters.q)$dominant_ctss - 1
	dominant_ctss_end = mcols(clusters.q)$dominant_ctss
	domAtStart <- mcols(clusters.q)$dominant_ctss == start(clusters.q)
	dominant_ctss_start[domAtStart] = start(clusters.q)[domAtStart]
	dominant_ctss_end[domAtStart]   = start(clusters.q)[domAtStart]
	
	sth = 0
	block_nr = rep(3, length(clusters.q))
	block_nr[which((dominant_ctss_end <= q.low.pos) | (dominant_ctss_end > q.up.pos))] = 4
	block_lengths = paste(1, formatC(q.up.pos - q.low.pos, format = 'f', digits = 0), 1, sep = ',')
	block_lengths[which(dominant_ctss_end <= q.low.pos)] = paste(1, 1, formatC(q.up.pos[which(dominant_ctss_end <= q.low.pos)] - q.low.pos[which(dominant_ctss_end <= q.low.pos)], format = 'f', digits = 0), 1, sep = ',')
	block_lengths[which(dominant_ctss_end > q.up.pos)] = paste(1, formatC(q.up.pos[which(dominant_ctss_end > q.up.pos)] - q.low.pos[which(dominant_ctss_end > q.up.pos)], format = 'f', digits = 0), 1, 1, sep = ',')
	block_rel_pos = paste(0, formatC(q.low.pos - start, format = 'f', digits = 0), end - start - 1, sep = ',')
	block_rel_pos[which(dominant_ctss_end <= q.low.pos)] = paste(0, formatC(dominant_ctss_start[which(dominant_ctss_end <= q.low.pos)] - start[which(dominant_ctss_end <= q.low.pos)], format = 'f', digits = 0), formatC(q.low.pos[which(dominant_ctss_end <= q.low.pos)] - start[which(dominant_ctss_end <= q.low.pos)], format = 'f', digits = 0), end[which(dominant_ctss_end <= q.low.pos)] - start[which(dominant_ctss_end <= q.low.pos)] - 1, sep = ',')
	block_rel_pos[which(dominant_ctss_end > q.up.pos)] = paste(0, formatC(q.low.pos[which(dominant_ctss_end > q.up.pos)] - start[which(dominant_ctss_end > q.up.pos)], format = 'f', digits = 0), formatC(dominant_ctss_start[which(dominant_ctss_end > q.up.pos)] - start[which(dominant_ctss_end > q.up.pos)], format = 'f', digits = 0), end[which(dominant_ctss_end > q.up.pos)] - start[which(dominant_ctss_end > q.up.pos)] - 1, sep = ',')
	
    removeFirstBlock <- unlist(lapply(strsplit(block_rel_pos, split = ",", fixed = T), function(x) {x[1] == x[2]}))
    removeLastBlock <- unlist(lapply(as.list(c(1:length(end))), function(x) {s <- tail(strsplit(block_rel_pos[x], split = ",", fixed = T)[[1]],2)[1]; l <- tail(strsplit(block_lengths[x], split = ",", fixed = T)[[1]],2)[1]; return((as.integer(end[x])-as.integer(start[x])) == (as.integer(s)+as.integer(l)))}))
    
    block_nr[removeFirstBlock] <- block_nr[removeFirstBlock] - 1
    block_nr[removeLastBlock] <- block_nr[removeLastBlock] - 1

    block_lengths[removeFirstBlock] <- sub("[[:digit:]]+,", "", block_lengths[removeFirstBlock])
    block_lengths[removeLastBlock] <- reverse(sub("[[:digit:]]+,", "", reverse(block_lengths[removeLastBlock])))
    
    block_rel_pos[removeFirstBlock] <- sub("[[:digit:]]+,", "", block_rel_pos[removeFirstBlock])
    block_rel_pos[removeLastBlock] <- reverse(sub("[[:digit:]]+,", "", reverse(block_rel_pos[removeLastBlock])))
    
		if(itemRgb){
            write(paste('track name="', track.name,'" description="', track.description,'" visibility="pack"', ' itemRgb="On"', sep = ''), file = track.file, append = app)
			write.table(data.frame(chr, formatC(start, format = 'f', digits = 0), formatC(end, format = 'f', digits = 0), name, score = rep(0,length(start)), strand, formatC(dominant_ctss_start,  format = 'f', digits = 0), formatC(dominant_ctss_end, format = 'f', digits = 0), cols, block_nr, block_lengths, block_rel_pos), file = track.file, append = T, col.names = F, row.names = F, quote = F, sep = '\t')
		}else{
            write(paste('track name="', track.name,'" description="', track.description,'" visibility="pack" color=', cols, sep = ''), file = track.file, append = app)
			write.table(data.frame(chr, formatC(start, format = 'f', digits = 0), formatC(end, format = 'f', digits = 0), name, score = rep(0,length(start)), strand, formatC(dominant_ctss_start,  format = 'f', digits = 0), formatC(dominant_ctss_end, format = 'f', digits = 0), rep(0, length(start)), block_nr, block_lengths, block_rel_pos), file = track.file, append = T, col.names = F, row.names = F, quote = F, sep = '\t')		
        }
		
	}else{

		if(itemRgb){
			write(paste('track name="', track.name,'" description="', track.description,'" visibility="pack"', ' itemRgb="On"', sep = ''), file = track.file, append = app)
			write.table(data.frame(chr, formatC(start, format = 'f', digits = 0), formatC(end, format = 'f', digits = 0), name, score = rep(0,length(start)), strand, formatC(start, format = 'f', digits = 0), formatC(end, format = 'f', digits = 0), cols), file = track.file, append = T, col.names = F, row.names = F, quote = F, sep = '\t')
		}else{
			write(paste('track name="', track.name,'" description="', track.description,'" visibility="pack" color=', cols, sep = ''), file = track.file, append = app)
			write.table(data.frame(chr, formatC(start, format = 'f', digits = 0), formatC(end, format = 'f', digits = 0), name, score = rep(0,length(start)), strand, formatC(start, format = 'f', digits = 0), formatC(end, format = 'f', digits = 0), rep(0, length(start))), file = track.file, append = T, col.names = F, row.names = F, quote = F, sep = '\t')		
		}
		
		
	}

}



