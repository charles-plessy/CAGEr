#' @include ClusteringMethods.R

###################################################################################################
# Implementation of an algorithm for correcting systematic G nucleotide addition bias to CAGE tags
# (as described in Carninci et al., Nature Genetics 2006, Supplementary Information, section 3-e)


.estimate.G.addition.and.correct <- function(ctss, G.chance, correction.orientation) {
	
	ctss.dt <- data.table(ctss)
	removedG <- pos <- NULL # Keep R CMD check happy
	ctss.count <- ctss.dt[, list(length(removedG), sum(removedG)), by = pos]
	setkeyv(ctss.count, cols = "pos")
	
	# select the 'G' positions in the genome that had some tags before moving the 'G' mismatch reads from previous position - these should be further corrected
	V1 <- V2 <- NULL # Keep R CMD check happy
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
	
	
	return(subset(ctss.final, ctss.final$nr_tags > 0))
	
}


#' .remove.added.G 
#' 
#' Non-exported private helper function
#' 
#' @noRd
#' @importFrom BSgenome getSeq

.remove.added.G <- function(reads.GRanges, genome, correctSystematicG = TRUE, sample.label) {

  reads.GRanges.plus  <- reads.GRanges[strand(reads.GRanges) == "+"]
  reads.GRanges.minus <- reads.GRanges[strand(reads.GRanges) == "-"]
	
	G.reads.plus <- which(substr(reads.GRanges.plus$seq, start = 1, stop = 1) == "G")
	G.reads.minus <- which(substr(reads.GRanges.minus$seq, start = reads.GRanges.minus$read.length, stop = reads.GRanges.minus$read.length) == "C")
	
	if(length(G.reads.plus)>0){
        G.mismatch.reads.plus <- G.reads.plus[getSeq(getRefGenome(genome), resize(reads.GRanges.plus[G.reads.plus], width = 1, fix = "start"), as.character = TRUE) != "G"]
        reads.GRanges.plus$removedG <- FALSE
        reads.GRanges.plus$removedG[G.mismatch.reads.plus] <- TRUE
        start(reads.GRanges.plus)[G.mismatch.reads.plus] <- start(reads.GRanges.plus)[G.mismatch.reads.plus] + as.integer(1)
        CTSS.plus <- data.frame(chr = as.character(seqnames(reads.GRanges.plus)), pos = start(reads.GRanges.plus), strand = "+", removedG = reads.GRanges.plus$removedG, stringsAsFactors = FALSE)
	}else{
		G.mismatch.reads.plus <- NULL
		CTSS.plus <- data.frame()
	}
	
	if(length(G.reads.minus)>0){
        G.mismatch.reads.minus <- G.reads.minus[getSeq(getRefGenome(genome), resize(reads.GRanges.minus[G.reads.minus], width = 1, fix = "start"), as.character = TRUE) != "G"]
        reads.GRanges.minus$removedG <- FALSE
        reads.GRanges.minus$removedG[G.mismatch.reads.minus] <- TRUE
        end(reads.GRanges.minus)[G.mismatch.reads.minus] <- end(reads.GRanges.minus)[G.mismatch.reads.minus] - as.integer(1)
        CTSS.minus <- data.frame(chr = as.character(seqnames(reads.GRanges.minus)), pos = end(reads.GRanges.minus), strand = "-", removedG = reads.GRanges.minus$removedG, stringsAsFactors = FALSE)
	}else{
		G.mismatch.reads.minus <- NULL
		CTSS.minus <- data.frame()
	}
	
	if(correctSystematicG){
		
		message("\t-> Estimating the frequency of adding a 'G' nucleotide and correcting the systematic bias...")
		
		# estimate chance of adding a 'G' nucleotide at the beginning of the CAGE tag
		# -> proportion of reads with unambigously added 'G' ('G' in the read and not in the genome) in all unambigous reads (reads with no 'G' at the beginning + reads with unambigously added 'G')
		G.chance <- (length(G.mismatch.reads.plus) + length(G.mismatch.reads.minus)) / ((length(reads.GRanges.plus) - length(G.reads.plus)) + (length(reads.GRanges.minus) - length(G.reads.minus)) + length(G.mismatch.reads.plus) + length(G.mismatch.reads.minus))
		
		f <- function(CTSS, filter, strand) {
      CTSS.G    <- CTSS[filter,]
    
      CTSS.G.corrected <-
        lapply( as.list(unique(CTSS.G$chr))
              , function(x) { ctss.corrected <-
                                .estimate.G.addition.and.correct(
                                     ctss                   = CTSS.G[CTSS.G$chr == x,]
                                   , G.chance               = G.chance
                                   , correction.orientation = ifelse(strand == "+", 1, -1))
                              ctss.corrected$chr <- x
                              return(ctss.corrected)
                             })
      CTSS.G.corrected <- do.call(rbind, CTSS.G.corrected)
      CTSS.G.corrected$strand <- strand
      CTSS.G.corrected <- CTSS.G.corrected[, c("chr", "pos", "strand", "nr_tags")]
      
      CTSS.no.G <- data.table(CTSS[-filter,])
      CTSS.no.G <- .byCtss(CTSS.no.G, "removedG", length)
      setnames(CTSS.no.G, c("chr", "pos", "strand", "nr_tags"))
      
      CTSS.final <- rbind(CTSS.G.corrected, as.data.frame(CTSS.no.G))
      CTSS.final <- .byCtss(data.table(CTSS.final), "nr_tags", sum)
      CTSS.final
    }
		
		if(nrow(CTSS.plus)>0){
			CTSS.plus.final <- f(CTSS.plus, G.reads.plus, "+")
		}else{
			CTSS.plus.final <- data.table()
		}

		if(nrow(CTSS.minus)>0){
			CTSS.minus.final <- f(CTSS.minus, G.reads.minus, "-")
		}else{
			CTSS.minus.final <- data.table()
		}
				
		CTSS <- data.table(rbind(as.data.frame(CTSS.plus.final), as.data.frame(CTSS.minus.final)))

	}else{
	
		CTSS <- rbind(CTSS.plus, CTSS.minus)
		CTSS <- CTSS[,c("chr", "pos", "strand")]
		CTSS$tag_count <- 1
		CTSS <- data.table(CTSS)
		CTSS <- .byCtss(CTSS, "tag_count", sum)
	}
  setnames(CTSS, c("chr", "pos", "strand", sample.label)) 
  setkeyv(CTSS, cols = c("chr", "pos", "strand"))

  gp <- CTSS(CTSS$chr, CTSS$pos, CTSS$strand, bsgenomeName = genome)
  score(gp) <- CTSS$score
  gp
}


#' .remove.added.G.CTSS
#' 
#' Non-exported private helper function
#' 
#' @examples 
#' gr <- GRanges("chr1", IRanges(1, 10), strand = c("+", "+", "-", "-"))
#' gr$seq <- c("ATTTAAATTT", "GTTTAAATTT", "TTTAAATTTA", "TTTAAATTTC")
#' gr$read.length <- 10
#' genome <- genomeName(exampleCAGEexp)
#' .remove.added.G.CTSS(gr, genome)
#' 
#' @noRd
#' @importFrom BSgenome getSeq

.remove.added.G.CTSS <- function(gr, genome, correctSystematicG = FALSE) {
    if(correctSystematicG == TRUE) stop ("correctSystematicG not supported in this function, use .remove.added.G instead.")
  gr <- promoters(gr, 0, 1)
  gr$genomeSeq <- getSeq(getRefGenome(genome), gr, as.character = TRUE)
  
  grl <- split(gr, strand(gr))

  removeOnPlus <- function(gr) {
    firstbase <- substr(gr$seq, start = 1, stop = 1)
    extraG <- firstbase == "G" & gr$genomeSeq != "G"
    ranges(gr[extraG]) <- shift(ranges(gr[extraG]), 1)
    gr
  }
  
  grl[["+"]] <- removeOnPlus(grl[["+"]])
  
  removeOnMinus <- function(gr) {
    firstbase <- substr(gr$seq,start = gr$read.length, stop = gr$read.length)
    extraG <- firstbase == "C" & gr$genomeSeq != "G"
    ranges(gr[extraG]) <- shift(ranges(gr[extraG]), -1)
    gr
  }
  
  grl[["-"]] <- removeOnMinus(grl[["-"]])

  CTSS(unlist(grl))
}
