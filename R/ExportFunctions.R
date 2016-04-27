###########################################################
# Functions for exporting results (graphical and textual)
#

.plotReverseCumulative <- function(values, col = "darkblue", title = "", col.title = "black", add = FALSE) {

	values <- values[values != 0]

# using data.table package
	v <- data.table(num = 1, nr_tags = values)
	v <- v[, sum(num), by = nr_tags]
	setkey(v, nr_tags)
	
#	v <- aggregate(values, by = list(values), FUN = length)
#	colnames(v) <- c('nr_tags', 'V1')
	
	if(add){
		lines(sort(c(v$nr_tags - (v$nr_tags - c(0, v$nr_tags[-length(v$nr_tags)]))/2, v$nr_tags + (c(v$nr_tags[-1], v$nr_tags[length(v$nr_tags)]+1) - v$nr_tags)/2)), rep(rev(cumsum(rev(v$V1))), each = 2), col = col, lwd = 2)
	}else{
		plot(sort(c(v$nr_tags - (v$nr_tags - c(0, v$nr_tags[-length(v$nr_tags)]))/2, v$nr_tags + (c(v$nr_tags[-1], v$nr_tags[length(v$nr_tags)]+1) - v$nr_tags)/2)), rep(rev(cumsum(rev(v$V1))), each = 2), xaxt = 'n', yaxt = 'n', log = "xy", type = "l", lwd = 2, col = col, xlab = "number of CAGE tags", ylab = "number of CTSSs (>= nr tags)", main = title, cex.axis = 1.8, cex.lab = 1.8, cex.main = 2.5, col.main = col.title, xlim = c(1, 10^5), ylim = c(1, 10^6))	
	
		axis(side = 1, at = as.integer(sapply(10^(seq(0,6,1)), function(x) {seq(x,10*x-1,x)})), labels = rep("", 63), cex.axis = 1.8)   		   
		axis(side = 1, at = 10^(seq(0,6,1)), labels = formatC(10^(seq(0,6,1)), format = "f", digits = 0), cex.axis = 1.8)		   
		axis(side = 2, at = as.integer(sapply(10^(seq(0,6,1)), function(x) {seq(x,10*x-1,x)})), labels = rep("", 63), cex.axis = 1.8)   		   
		axis(side = 2, at = 10^(seq(0,6,1)), labels = formatC(10^(seq(0,6,1)), format = "f", digits = 0), cex.axis = 1.8)
	}
	
}

.export.bedgraph <- function(data.rd, name, description, file_name, append = F) {
	data.ucsc <- as(data.rd, "UCSCData")
	data.ucsc@trackLine <- new("BasicTrackLine", name = name, description = description, visibility="full")	
	export(data.ucsc, con = file_name, format = "ucsc", subformat = "bedGraph", append = append)
}


.export.bw.all <- function(data, sample.labels, v, genome) {
    data.plus <- subset(data, strand == "+")
    data.minus <- subset(data, strand == "-")
    
    rd.list <- lapply(as.list(sample.labels), function(x) {
        lapply(as.list(c("plus", "minus")), function(y) {
            d <- get(paste("data.", y, sep = ""))
            d <- d[, c("chr", "pos", x)]
            colnames(d) <- c("chr", "pos", "score")
            if(nrow(d)>0){
                d <- subset(d, score>0)
                d.rd <- GRanges(seqnames = d$chr, ranges=IRanges(start = d$pos, end = d$pos), strand = "+", score = d$score, seqlengths=seqlengths(genome))
            }else{
                d.rd <- GRanges()
                seqlengths(d.rd) <- seqlengths(genome)
            }
            if(y == "minus"){
                d.rd$score <- -1 * d.rd$score
            }
            return(d.rd)
							 }
							 )
    }
    )
    names(rd.list) <- sample.labels

    a <- lapply(sample.labels, function(x){
            if(length(rd.list[[x]][[1]]) > 0){
                export.bw(rd.list[[x]][[1]], con = paste(x, ".CTSS.", v, ".plus.bw", sep = ""))
            }
            if(length(rd.list[[x]][[2]]) > 0){
                export.bw(rd.list[[x]][[2]], con = paste(x, ".CTSS.", v, ".minus.bw", sep = ""))
            }
        })
    
    description.lines = data.frame(description = unlist(lapply(sample.labels, function(x) {paste('track type=bigWig name="', paste(x, "CTSS", v, c('plus"', 'minus"'), sep = " "), ' description="', paste(x, " CTSS ", v, ' (', c("plus", "minus"), ' strand)"', sep = ""), " bigDataUrl=", paste(x, ".CTSS.", v, c(".plus", ".minus"), ".bw", sep = ""), sep = "")})))
    write.table(description.lines, file = paste("CTSS.", v, ".all.samples.track.description.txt", sep = ""), col.names = F, row.names = F, sep = "\t", quote = F)
}

.export.bedgraph.all <- function(data, sample.labels, v, oneFile) {
	data.plus <- subset(data, strand == "+")
	data.minus <- subset(data, strand == "-")
		
	rd.list <- lapply(as.list(sample.labels), function(x) {
					  lapply(as.list(c("plus", "minus")), function(y) {
							 d <- get(paste("data.", y, sep = ""))
							 d <- d[, c("chr", "pos", x)] 
							 colnames(d) <- c("chr", "pos", "score")
							 if(nrow(d)>0){
							 d <- subset(d, score>0)
							 d.rd <- GRanges(seqnames = d$chr, ranges=IRanges(start = d$pos, end = d$pos), strand = "+", score = d$score)
							 }else{
							 d.rd <- GRanges()
							 }
							 if(y == "minus"){
							 d.rd$score <- -1 * d.rd$score
							 }
							 return(d.rd)
							 }
							 )
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

.plot.clusters.beanplots <- function(value.matrix, cl, cl.method, dim.som.x, dim.som.y, ylim = c(0,2), las = 0, labels = colnames(value.matrix), titles = 'number', cex.axis = 1, cex.main = 5, cex.lab = 3, cols = c("red", "gold", "green", "blue")) {
	
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



.make.cluster.bed.track <- function(clusters.q, use.blocks = T, q.low = NULL, q.up = NULL, track.file, track.name = 'q_track', track.description = 'q track', name = ".", cols = "0,0,0", itemRgb = FALSE, app = T) {
	
	if(!(app) & file.exists(track.file)){
		file.remove(track.file)
	}	

	chr = clusters.q$chr
	start = clusters.q$start
	end = clusters.q$end
	strand = clusters.q$strand

	if(use.blocks){
		
	q.low.pos = clusters.q[,paste('qLow_', q.low, sep = '')] - 1 
	q.up.pos = clusters.q[,paste('qUp_', q.up, sep = '')]
		
	dominant_ctss_start = clusters.q$dominant_ctss - 1
	dominant_ctss_end = clusters.q$dominant_ctss
	dominant_ctss_start[clusters.q$dominant_ctss == clusters.q$start] = clusters.q$start[clusters.q$dominant_ctss == clusters.q$start]
	dominant_ctss_end[clusters.q$dominant_ctss == clusters.q$start] = clusters.q$start[clusters.q$dominant_ctss == clusters.q$start]
	
	sth = 0
	block_nr = rep(3, nrow(clusters.q))
	block_nr[which((clusters.q$dominant_ctss <= q.low.pos) | (clusters.q$dominant_ctss > q.up.pos))] = 4
	block_lengths = paste(1, formatC(q.up.pos - q.low.pos, format = 'f', digits = 0), 1, sep = ',')
	block_lengths[which(clusters.q$dominant_ctss <= q.low.pos)] = paste(1, 1, formatC(q.up.pos[which(clusters.q$dominant_ctss <= q.low.pos)] - q.low.pos[which(clusters.q$dominant_ctss <= q.low.pos)], format = 'f', digits = 0), 1, sep = ',')
	block_lengths[which(clusters.q$dominant_ctss > q.up.pos)] = paste(1, formatC(q.up.pos[which(clusters.q$dominant_ctss > q.up.pos)] - q.low.pos[which(clusters.q$dominant_ctss > q.up.pos)], format = 'f', digits = 0), 1, 1, sep = ',')
	block_rel_pos = paste(0, formatC(q.low.pos - start, format = 'f', digits = 0), end - start - 1, sep = ',')
	block_rel_pos[which(clusters.q$dominant_ctss <= q.low.pos)] = paste(0, formatC(dominant_ctss_start[which(clusters.q$dominant_ctss <= q.low.pos)] - start[which(clusters.q$dominant_ctss <= q.low.pos)], format = 'f', digits = 0), formatC(q.low.pos[which(clusters.q$dominant_ctss <= q.low.pos)] - start[which(clusters.q$dominant_ctss <= q.low.pos)], format = 'f', digits = 0), end[which(clusters.q$dominant_ctss <= q.low.pos)] - start[which(clusters.q$dominant_ctss <= q.low.pos)] - 1, sep = ',')
	block_rel_pos[which(clusters.q$dominant_ctss > q.up.pos)] = paste(0, formatC(q.low.pos[which(clusters.q$dominant_ctss > q.up.pos)] - start[which(clusters.q$dominant_ctss > q.up.pos)], format = 'f', digits = 0), formatC(dominant_ctss_start[which(clusters.q$dominant_ctss > q.up.pos)] - start[which(clusters.q$dominant_ctss > q.up.pos)], format = 'f', digits = 0), end[which(clusters.q$dominant_ctss > q.up.pos)] - start[which(clusters.q$dominant_ctss > q.up.pos)] - 1, sep = ',')
	
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



