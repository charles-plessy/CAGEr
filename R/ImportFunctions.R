###################################################################################################
# Implementation of an algorithm for correcting systematic G nucleotide addition bias to CAGE tags
# (as described in Carninci et al., Nature Genetics 2006, Supplementary Information, section 3-e)


.estimate.G.addition.and.correct <- function(ctss, G.chance, correction.orientation) {
	
	ctss.dt <- data.table(ctss)
	ctss.count <- ctss.dt[, list(length(removedG), sum(removedG)), by = pos]
	setkey(ctss.count, pos)
	
	# select the 'G' positions in the genome that had some tags before moving the 'G' mismatch reads from previous position - these should be further corrected
	ctss.to.correct <- data.frame(subset(ctss.count, V1 != V2))
	
	if(nrow(ctss.to.correct) > 0){
	# iterate sequentially through 'G' positions and correct for the number of reads that have to be moved to downstream position
	
	if(correction.orientation > 0){
		ctss.gap <- c(Inf, diff(ctss.to.correct$pos))
	}else if(correction.orientation < 0){
		ctss.gap <- c(diff(ctss.to.correct$pos), Inf)
	}
	G.start <- which(ctss.gap != 1)
	G.follow <- which(ctss.gap == 1)
	ctss.to.append <- data.frame()
	
	while(length(G.start) > 0) {
		
		F <- as.integer(pmax(round(ctss.to.correct$V1[G.start] - ctss.to.correct$V2[G.start]/G.chance), 0))
		ctss.to.correct$V1[G.start] <- ctss.to.correct$V1[G.start] - F
		idx <- G.start + correction.orientation
		if(correction.orientation > 0){
			idx[idx == (nrow(ctss.to.correct) + 1)] <- 1			
		}else if(correction.orientation < 0){
			idx[idx == 0] <- nrow(ctss.to.correct)
		}
		G.start.followed <- ctss.to.correct$pos[idx] %in% ctss.to.correct$pos[G.follow]
		ctss.to.correct$V1[G.start[G.start.followed] + correction.orientation] <- ctss.to.correct$V1[G.start[G.start.followed] + correction.orientation] + F[G.start.followed]
		ctss.to.correct$V2[G.start[G.start.followed] + correction.orientation] <- F[G.start.followed]
		ctss.to.append <- rbind(ctss.to.append, data.frame(pos = ctss.to.correct$pos[G.start[!G.start.followed]] + correction.orientation, V1 = F[!G.start.followed], V2 = F[!G.start.followed]))
	
		G.start <- (G.start + correction.orientation)[G.start.followed]
	
	}
	
	ctss.final <- rbind(ctss.to.correct, ctss.to.append)
	ctss.final <- rbind(ctss.final, as.data.frame(ctss.count[V1 == V2]))
	ctss.final <- ctss.final[order(ctss.final$pos),]
	ctss.final <- data.frame(pos = ctss.final$pos, nr_tags = ctss.final$V1)
	
	}else{
		ctss.final <- data.frame(pos = ctss.count$pos, nr_tags = ctss.count$V1)
	}
	
	return(subset(ctss.final, nr_tags > 0))
	
}


.remove.added.G <- function(reads.GRanges.plus, reads.GRanges.minus, genome, correctSystematicG = TRUE) {

	message("\t-> Removing the first base of the reads if 'G' and not aligned to the genome...")
	
	G.reads.plus <- which(substr(elementMetadata(reads.GRanges.plus)$seq, start = 1, stop = 1) == "G")
	G.reads.minus <- which(substr(elementMetadata(reads.GRanges.minus)$seq, start = elementMetadata(reads.GRanges.minus)$read.length, stop = elementMetadata(reads.GRanges.minus)$read.length) == "C")
	
	if(length(G.reads.plus)>0){
	G.mismatch.reads.plus <- G.reads.plus[getSeq(genome, resize(reads.GRanges.plus[G.reads.plus], width = 1, fix = "start"), as.character = TRUE) != "G"]
	elementMetadata(reads.GRanges.plus)$removedG <- FALSE
	elementMetadata(reads.GRanges.plus)$removedG[G.mismatch.reads.plus] <- TRUE
	start(reads.GRanges.plus)[G.mismatch.reads.plus] <- start(reads.GRanges.plus)[G.mismatch.reads.plus] + as.integer(1)
	CTSS.plus <- data.frame(chr = as.character(seqnames(reads.GRanges.plus)), pos = start(reads.GRanges.plus), strand = "+", removedG = elementMetadata(reads.GRanges.plus)$removedG, stringsAsFactors = FALSE)
	}else{
		G.mismatch.reads.plus <- NULL
		CTSS.plus <- data.frame()
	}
	
	if(length(G.reads.minus)>0){
	G.mismatch.reads.minus <- G.reads.minus[getSeq(genome, resize(reads.GRanges.minus[G.reads.minus], width = 1, fix = "start"), as.character = TRUE) != "G"]
	elementMetadata(reads.GRanges.minus)$removedG <- FALSE
	elementMetadata(reads.GRanges.minus)$removedG[G.mismatch.reads.minus] <- TRUE
	end(reads.GRanges.minus)[G.mismatch.reads.minus] <- end(reads.GRanges.minus)[G.mismatch.reads.minus] - as.integer(1)
	CTSS.minus <- data.frame(chr = as.character(seqnames(reads.GRanges.minus)), pos = end(reads.GRanges.minus), strand = "-", removedG = elementMetadata(reads.GRanges.minus)$removedG, stringsAsFactors = FALSE)
	}else{
		G.mismatch.reads.minus <- NULL
		CTSS.minus <- data.frame()
	}
	
	if(correctSystematicG){
		
		message("\t-> Estimating the frequency of adding a 'G' nucleotide and correcting the systematic bias...")
		
		# estimate chance of adding a 'G' nucleotide at the beginning of the CAGE tag
		# -> proportion of reads with unambigously added 'G' ('G' in the read and not in the genome) in all unambigous reads (reads with no 'G' at the beginning + reads with unambigously added 'G')
		G.chance <- (length(G.mismatch.reads.plus) + length(G.mismatch.reads.minus)) / ((length(reads.GRanges.plus) - length(G.reads.plus)) + (length(reads.GRanges.minus) - length(G.reads.minus)) + length(G.mismatch.reads.plus) + length(G.mismatch.reads.minus))
		
		if(nrow(CTSS.plus)>0){

			CTSS.G.plus <- CTSS.plus[G.reads.plus,]
			CTSS.G.plus.corrected <- lapply(as.list(unique(CTSS.G.plus$chr)), function(x) {ctss.corrected <- .estimate.G.addition.and.correct(ctss = subset(CTSS.G.plus, chr == x), G.chance = G.chance, correction.orientation = 1); ctss.corrected$chr = x; return(ctss.corrected)})
			CTSS.G.plus.corrected <- do.call(rbind, CTSS.G.plus.corrected)
			CTSS.G.plus.corrected$strand <- "+"
			CTSS.G.plus.corrected <- CTSS.G.plus.corrected[,c("chr", "pos", "strand", "nr_tags")]

			CTSS.no.G.plus <- data.table(CTSS.plus[-G.reads.plus,])
			CTSS.no.G.plus <- CTSS.no.G.plus[, length(removedG), by = list(chr, pos, strand)]
			setnames(CTSS.no.G.plus, c("chr", "pos", "strand", "nr_tags"))
			CTSS.plus.final <- rbind(CTSS.G.plus.corrected, as.data.frame(CTSS.no.G.plus))
			CTSS.plus.final <- data.table(CTSS.plus.final)
			CTSS.plus.final <- CTSS.plus.final[, sum(nr_tags), by = list(chr, pos, strand)]

		}else{
		
			CTSS.plus.final <- data.table()
			
		}

		if(nrow(CTSS.minus)>0){

			CTSS.G.minus <- CTSS.minus[G.reads.minus,]
			CTSS.G.minus.corrected <- lapply(as.list(unique(CTSS.G.minus$chr)), function(x) {ctss.corrected <- .estimate.G.addition.and.correct(ctss = subset(CTSS.G.minus, chr == x), G.chance = G.chance, correction.orientation = -1); ctss.corrected$chr = x; return(ctss.corrected)})
			CTSS.G.minus.corrected <- do.call(rbind, CTSS.G.minus.corrected)
			CTSS.G.minus.corrected$strand <- "-"
			CTSS.G.minus.corrected <- CTSS.G.minus.corrected[,c("chr", "pos", "strand", "nr_tags")]
		
			CTSS.no.G.minus <- data.table(CTSS.minus[-G.reads.minus,])
			CTSS.no.G.minus <- CTSS.no.G.minus[, length(removedG), by = list(chr, pos, strand)]
			setnames(CTSS.no.G.minus, c("chr", "pos", "strand", "nr_tags"))
			CTSS.minus.final <- rbind(CTSS.G.minus.corrected, as.data.frame(CTSS.no.G.minus))
			CTSS.minus.final <- data.table(CTSS.minus.final)
			CTSS.minus.final <- CTSS.minus.final[, sum(nr_tags), by = list(chr, pos, strand)]

		}else{
			
			CTSS.minus.final <- data.table()
			
		}
				
		CTSS <- data.table(rbind(as.data.frame(CTSS.plus.final), as.data.frame(CTSS.minus.final)))

	}else{
	
		CTSS <- rbind(CTSS.plus, CTSS.minus)
		CTSS <- CTSS[,c("chr", "pos", "strand")]
		CTSS$tag_count <- 1
		CTSS <- data.table(CTSS)
		CTSS <- CTSS[, as.integer(sum(tag_count)), by = list(chr, pos, strand)]		
	
	}
	
	return(CTSS)
	
}


