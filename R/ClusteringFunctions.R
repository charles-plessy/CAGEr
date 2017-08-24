#' @include CTSS.R

####################################################################################
# Implementation of simple distance-based clustering of data attached to sequences
#

#' @name distclu-functions
#' 
#' @title Private functions for distance clustering.
#' 
#' @description The flow of data is that a \code{GRanges} object of CTSSes is progressively
#' deconstructed, and data to form the clusters is progressively integrated in a
#' \code{\link{data.table}} object, which is finally converted to GRanges at the end.  Doing
#' the whole clustering with GRanges is more elegant, but looping on a GRangesList was just
#' too slow.  Maybe the operation on the \code{data.table} is more efficient because it is
#' vectorised.
NULL

#' @name .cluster.ctss.strand
#' @rdname distclu-functions
#' 
#' @description \code{.cluster.ctss.strand} does the strandless distance clustering of strandless
#' CTSS positions from a single chromosome.  Input does not need to be sorted, but \emph{pay
#' attention that the output is sorted}.
#' 
#' @param ctss.iranges.chr A IRanges object.
#' @param max.dist See \code{\link{clusterCTSS}}.
#' 
#' @return \code{.cluster.ctss.strand} returns an \code{\link{data.table}} object containing
#'  arbitrary cluster IDs (as integers) for each CTSS.
#'  
#' @importFrom data.table data.table
#' 
#' @examples
#' 
#' #.cluster.ctss.strand
#' ctss.iranges.chr <- IRanges(c(1,3,4,12,14,25,28), w=1)
#' .cluster.ctss.strand(ctss.iranges.chr, 5)
#' 
#' ce <- readRDS(system.file(package = "CAGEr", "extdata/CAGEexp.rds"))
#' ctss.chr <- CTSScoordinatesGR(ce)
#' ctss.chr <- ctss.chr[strand(ctss.chr) == "+"]
#' ctss.iranges.chr <- ranges(ctss.chr)
#' # Same result if not sorted
#' identical(
#'   .cluster.ctss.strand(ctss.iranges.chr, 20),
#'   .cluster.ctss.strand(ctss.iranges.chr[sample(seq_along(ctss.iranges.chr))], 20)
#' )
#' # Returns an emtpy data.table object if given an empty IRanges object.
#' identical(.cluster.ctss.strand(IRanges(), 20), data.table())

setGeneric(".cluster.ctss.strand", function(ctss.iranges.chr, max.dist) standardGeneric(".cluster.ctss.strand"))

setMethod(".cluster.ctss.strand", "IRanges", function(ctss.iranges.chr, max.dist) {
  if (identical(ctss.iranges.chr, IRanges())) return(data.table())
  v <- Rle(rep(TRUE, max(start(ctss.iranges.chr)) + max.dist))
  v[start(ctss.iranges.chr) + max.dist] <- FALSE
  longRuns <- which(runLength(v) >= max.dist & runValue(v))
  c.starts <- cumsum(runLength(v))[longRuns] - max.dist
  c.ends <- c(cumsum(runLength(v))[(longRuns - 1)[-1]], length(v)) - max.dist
  clusters <- IRanges(start = c.starts, end = c.ends)
  o <- findOverlaps(clusters, ctss.iranges.chr)
  data.table(id = queryHits(o))
})


#' @name .cluster.ctss.chr
#' @rdname distclu-functions
#' 
#' @description \code{.cluster.ctss.chr} does the stranded distance clustering of CTSS on a
#' single chromosome, by dispatching both strands to \code{.cluster.ctss.strand} and merging
#' the results, taking care keep the cluster IDs unique.  Be careful that this function does
#' not look at the score.
#' 
#' @param ctss.chr A CTSS.chr object.
#' @param max.dist See \code{\link{clusterCTSS}}.
#' 
#' @return \code{.cluster.ctss.chr} returns a \code{\link{data.table}} object representing the
#' chromosome coordinates (\code{chr}, \code{pos}, \code{strand}) of each CTSS, with their
#' cluster ID (\code{id}).
#' 
#' @importFrom data.table data.table
#' @importFrom S4Vectors decode
#' 
#' @examples 
#' 
#' #.cluster.ctss.chr
#' ce <- readRDS(system.file(package = "CAGEr", "extdata/CAGEexp.rds"))
#' ctss.chr <- CAGEr:::.CTSS.chr(CTSScoordinatesGR(ce))
#' .cluster.ctss.chr(ctss.chr, 20)

setGeneric(".cluster.ctss.chr", function(ctss.chr, max.dist) standardGeneric(".cluster.ctss.chr"))

setMethod(".cluster.ctss.chr", "CTSS.chr", function(ctss.chr, max.dist) {
  clusters.p <- .cluster.ctss.strand(ranges(ctss.chr[strand(ctss.chr) == "+"]), max.dist)
  clusters.m <- .cluster.ctss.strand(ranges(ctss.chr[strand(ctss.chr) == "-"]), max.dist)
  clusters.m <- clusters.m + max(0, suppressWarnings(max(clusters.p)))
  data.table( rbind(clusters.p, clusters.m)
            , chr    = as.character(seqnames(ctss.chr))
            , pos    = start(ctss.chr)
            , strand = decode(strand(ctss.chr)))
})


#' @name .ctss2clusters
#' @rdname distclu-functions
#' 
#' @description \code{.ctss2clusters} does the stranded distance clustering of CTSS.
#' 
#' @param ctss A CTSS object with a score column.
#' @param max.dist See \code{\link{clusterCTSS}}.
#' @param useMulticore,nrCores See clusterCTSS.
#' 
#' @return \code{.ctss2clusters} returns a \code{\link{data.table}} object representing the
#' cluster ID (\code{id}), chromosome coordinates (\code{chr}, \code{pos}, \code{strand}) and
#' the \code{score} of each CTSS.
#' 
#' @importFrom data.table data.table
#' @importFrom S4Vectors decode
#' 
#' @examples 
#' 
#' # .ctss2clusters
#' ce <- readRDS(system.file(package = "CAGEr", "extdata/CAGEexp.rds"))
#' normalizeTagCount(ce)
#' ctss <- CAGEr:::.CTSS(CTSScoordinatesGR(ce))
#' score(ctss) <- CTSSnormalizedTpmDF(ce)[[1]]
#' seqnames(ctss)[rep(c(T,F), length(ctss) / 2)] <- "chr16"
#' ctss
#' .ctss2clusters(ctss, 20)

setGeneric(".ctss2clusters", function(ctss, max.dist = 20, useMulticore = FALSE, nrCores = NULL) standardGeneric(".ctss2clusters"))

setMethod(".ctss2clusters", "CTSS", function(ctss, max.dist, useMulticore, nrCores) {
  ctss <- sort(ctss)
  ctss <- ctss[score(ctss) != 0]
  ctss.list <- split(ctss, droplevels(seqnames(ctss)))
  ctss.list <- lapply(ctss.list, CAGEr:::.CTSS.chr)
  
  if(useMulticore == TRUE){
    requireNamespace("parallel")
    if(is.null(nrCores)){
      nrCores <- detectCores()
    }
		ctss.list <- mclapply(ctss.list, .cluster.ctss.chr, max.dist = max.dist, mc.cores = nrCores)		
	}else{
		ctss.list <- lapply(ctss.list, .cluster.ctss.chr, max.dist = max.dist)
	}
	max.clust <- sapply(ctss.list, function(x) {max(x$id)})
	max.clust <- cumsum(c(0, max.clust[-length(max.clust)]))
	
  for (x in seq_along(max.clust))
    ctss.list[[x]]$id <- ctss.list[[x]]$id + max.clust[x]
	
	dt <- do.call(rbind, ctss.list)
	dt$tpm <- decode(score(ctss))
	dt
})


#' @name .summarize.clusters
#' @rdname distclu-functions
#' @description \code{.summarize.clusters} calculates the number of CTSS, and
#' the position and score of a main peak, for each cluster.
#' 
#' @param ctss.clustered A \code{\link{data.table}} object representing the cluster ID
#'        (\code{id}), chromosome coordinates (\code{chr}, \code{pos}, \code{strand}) and
#'        the \code{score} of each CTSS.
#' @param removeSingletons Remove \dQuote{singleton} clusters that span only a single
#'        nucleotide ? (default = FALSE).
#' @param keepSingletonsAbove Even if \code{removeSingletons = TRUE}, keep singletons when
#'        their score is aboove threshold (default = \code{Inf}).
#' @return \code{.summarize.clusters} returns GRanges describing the clusters.
#' 
#' @importFrom data.table setnames
#' 
#' @examples 
#' 
#' # .summarize.clusters
#' ce <- readRDS(system.file(package = "CAGEr", "extdata/CAGEexp.rds"))
#' normalizeTagCount(ce)
#' ctss <- CAGEr:::.CTSS(CTSScoordinatesGR(ce))
#' score(ctss) <- CTSSnormalizedTpmDF(ce)[[1]]
#' seqnames(ctss)[rep(c(T,F), length(ctss) / 2)] <- "chr16"
#' clusters <- .ctss2clusters(ctss, 20)
#' .summarize.clusters(clusters)
#' .summarize.clusters(clusters, removeSingletons = TRUE)
#' .summarize.clusters(clusters, removeSingletons = TRUE, keepSingletonsAbove = 5)

setGeneric(".summarize.clusters", function(ctss.clustered, max.dist = 20, removeSingletons = FALSE, keepSingletonsAbove = Inf) standardGeneric(".summarize.clusters"))

setMethod(".summarize.clusters", "data.table", function(ctss.clustered, max.dist, removeSingletons, keepSingletonsAbove) {

  find.dominant.idx <- function (x) {
    w <- which(x == max(x))
    w[ceiling(length(w)/2)]
  }
	
  clusters <- ctss.clustered[ , list( chr[1]
                                    , min(pos)
                                    , max(pos)
                                    , strand[1]
                                    , length(pos)
                                    , pos[find.dominant.idx(tpm)]
                                    , sum(tpm)
                                    , max(tpm))
                              , by = id]
  setnames(clusters, c( "cluster"
                      , "chr", "start", "end", "strand"
                      , "nr_ctss", "dominant_ctss", "tpm", "tpm.dominant_ctss"))
  
	if(removeSingletons)
		clusters <- subset(clusters, nr_ctss > 1 | tpm >= keepSingletonsAbove)
  
  gr <- GRanges( seqnames = Rle(factor(clusters$chr))
               , ranges   = IRanges(clusters$start, clusters$end)
               , strand   = clusters$strand
               , score    = Rle(clusters$tpm)
               , nr_ctss  = clusters$nr_ctss
               , dominant_ctss = clusters$dominant_ctss
               , tpm.dominant_ctss = Rle(clusters$tpm.dominant_ctss)
  )
  names(gr) <- seq_along(gr)
  gr
})


#' @name .distclu
#' @rdname distclu-functions
#' @description  \code{.distclu} receives the data from the main \code{clusterCTSS} and
#' dispatches each for (possibly parallel) processing.
#' 
#' @param se A \code{\link{SummarizedExperiment}} object representing the CTSSes and
#'        their expression in each sample.
#' @param max.dist See \code{\link{clusterCTSS}}.
#' @param removeSingletons Remove \dQuote{singleton} clusters that span only a single
#'        nucleotide ? (default = FALSE).
#' @param keepSingletonsAbove Even if \code{removeSingletons = TRUE}, keep singletons when
#'        their score is aboove threshold (default = \code{Inf}).
#' @param useMulticore,nrCores See clusterCTSS.
#'        
#' @return \code{.distclu} returns GRanges describing the clusters.
#' 
#' @importFrom SummarizedExperiment rowRanges
#' 
#' @examples 
#' 
#' # .distclu
#' ce <- readRDS(system.file(package = "CAGEr", "extdata/CAGEexp.rds"))
#' normalizeTagCount(ce)
#' .distclu(CTSStagCountSE(ce))
#' \dontrun{
#' .distclu(CTSStagCountSE(ce), useMulticore = TRUE)
#' }
#' 
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' .distclu(CTSStagCountSE(exampleCAGEset))

setGeneric(".distclu", function(se, max.dist = 20, removeSingletons = FALSE, keepSingletonsAbove = Inf, useMulticore = FALSE, nrCores = NULL) standardGeneric(".distclu"))

setMethod(".distclu", "SummarizedExperiment", function(se, max.dist, removeSingletons, keepSingletonsAbove, useMulticore, nrCores) {
	
  ctss.cluster.list <- list()
  for(s in colnames(se)) {
    message("\t-> ", s)
    d <- CAGEr:::.CTSS(rowRanges(se))
    score(d) <- assays(se)[["normalizedTpmMatrix"]][[s]]
    d <- subset(d, score(d) > 0)
    clusters <- .ctss2clusters(ctss = d, max.dist = max.dist, useMulticore = useMulticore, nrCores = nrCores)
    ctss.cluster.list[[s]] <- clusters
  }

  if (useMulticore == TRUE) {
    requireNamespace("parallel")
    if(is.null(nrCores)){
      nrCores <- detectCores()
    }
    ctss.cluster.list <- mclapply(ctss.cluster.list, .summarize.clusters, removeSingletons = removeSingletons, keepSingletonsAbove = keepSingletonsAbove, mc.cores = nrCores)
  } else {
    ctss.cluster.list <- lapply(ctss.cluster.list, .summarize.clusters, removeSingletons = removeSingletons, keepSingletonsAbove = keepSingletonsAbove)
}
  GRangesList(ctss.cluster.list)
})


#################################################################################################################
# Implementation of Paraclu - parametric clustering of data attached to sequences (http://www.cbrc.jp/paraclu/) 
# Reference: A code for transcription initiation in mammalian genomes, 
# MC Frith, E Valen, A Krogh, Y Hayashizaki, P Carninci, A Sandelin, Genome Research 2008 18(1):1-12)
#

.paraclu1 <- function(ctss) {
	
	sit <- nrow(ctss)
	tot <- sum(ctss$tpm)
	if(sit == 1) {
		min_density <- Inf
		br <- NA
	}else{
		densities_forward <- cumsum(ctss$tpm)[-sit]/(ctss$pos[2:sit] - ctss$pos[1])
		densities_reverse <- cumsum(rev(ctss$tpm))[-sit]/(ctss$pos[sit] - ctss$pos[(sit-1):1])
		min_densities = c(min(densities_forward), min(densities_reverse))
		breaks <- c(which(densities_forward == min_densities[1])[1] + 1, sit + 1 - which(densities_reverse == min_densities[2])[1])
		min_density <- min(min_densities)
		br <- breaks[tail(which(min_densities == min_density),1)]
	}
	
	return(as.list(c(br, min_density, tot, sit)))
}


.paraclu2 <- function(ctss, min_density = -Inf, clusters.df = data.frame()) {
	
	if(nrow(ctss)>0){
		
		ctss <- ctss[order(ctss$pos),]
		params <- .paraclu1(ctss)
		br <- params[[1]]
		max_density <- params[[2]]
		tot <- params[[3]]
		sit <- params[[4]]
		
		if(!(max_density == Inf)){
			new_min <- max(min_density, max_density)
			clusters.df <- rbind(.paraclu2(ctss = ctss[1:(br-1),], min_density = new_min, clusters.df = clusters.df), .paraclu2(ctss = ctss[br:nrow(ctss),], min_density = new_min, clusters.df = clusters.df))
		}
		
		return(rbind(clusters.df, data.frame(chr = ctss$chr[1], start = min(ctss$pos), end = max(ctss$pos), strand = ctss$strand[1], nr_ctss = sit, dominant_ctss = ctss$pos[which(ctss$tpm == max(ctss$tpm))[ceiling(length(which(ctss$tpm == max(ctss$tpm)))/2)]], tpm = tot, tpm.dominant_ctss = ctss$tpm[which(ctss$tpm == max(ctss$tpm))[ceiling(length(which(ctss$tpm == max(ctss$tpm)))/2)]], min_d = min_density, max_d= max_density)))
		
	}else{
		return(clusters.df)
	}
	
}


.paraclu3 <- function(ctss.df, minStability = 1, maxLength = 500, removeSingletons = FALSE, keepSingletonsAbove = Inf, reduceToNonoverlapping = TRUE, useMulticore = FALSE, nrCores = NULL){
	
	if(useMulticore == TRUE){
		library(parallel)
		if(is.null(nrCores)){
			nrCores <- detectCores()
		}		
		
		ctss.df.plus.list <-lapply(as.list(unique(ctss.df$chr)), function(x) {subset(ctss.df, chr == x & strand == "+")})
		ctss.df.minus.list <-lapply(as.list(unique(ctss.df$chr)), function(x) {subset(ctss.df, chr == x & strand == "-")})
		ctss.df.list <- append(ctss.df.plus.list, ctss.df.minus.list)
		clusters.list <- mclapply(ctss.df.list, .paraclu2, mc.cores = nrCores)
		n <- length(clusters.list)/2
		clusters.list <- lapply(as.list(c(1:n)), function(x) {rbind(clusters.list[[x]], clusters.list[[x+n]])})
	
	}else{
	
		clusters.list <- list()
		for(ch in unique(ctss.df$chr)){
			cl.plus <- .paraclu2(ctss = subset(ctss.df, chr == ch & strand == "+"))
			cl.minus <- .paraclu2(ctss = subset(ctss.df, chr == ch & strand == "-"))
			clusters.list[[ch]] <- rbind(cl.plus, cl.minus)
		}
	}
	
	clusters <- do.call(rbind, clusters.list)
	
	clusters <- subset(clusters, (max_d >= (minStability * min_d)) & ((end - start + 1) <= maxLength))
	
	if(removeSingletons == TRUE){
		clusters <- subset(clusters, !(start == end & tpm < keepSingletonsAbove))
	}
	
	if(reduceToNonoverlapping == TRUE){
		clusters.gr <- GRanges(seqnames = clusters$chr, ranges = IRanges(start = clusters$start, end = clusters$end), strand = clusters$strand, elementMetadata = clusters)
		o <- findOverlaps(clusters.gr, drop.self = TRUE, type = "within")
		clusters.gr <- clusters.gr[-queryHits(o)]
		clusters <- subset(clusters, paste(chr, strand, start, end, sep = ".") %in% paste(seqnames(clusters.gr), strand(clusters.gr), start(clusters.gr), end(clusters.gr), sep = "."))
	}
	clusters <- cbind(cluster = c(1:nrow(clusters)), clusters)
	clusters$start <- clusters$start - 1
	rownames(clusters) <- c(1:nrow(clusters))
	colnames(clusters)[which(colnames(clusters) == "min_d")] <- "min_density"
	colnames(clusters)[which(colnames(clusters) == "max_d")] <- "max_density"
	
	return(clusters)
	
}


.paraclu <- function(data, sample.labels, minStability = 1, maxLength = 500, removeSingletons = FALSE, keepSingletonsAbove = Inf, reduceToNonoverlapping = TRUE, useMulticore = FALSE, nrCores = NULL) {
	
	ctss.cluster.list <- list()
	for(s in sample.labels) {
		
		message("\t-> ", s)
		d <- data[,c("chr", "pos", "strand", s)]
		colnames(d) <- c("chr", "pos", "strand", "tpm")
		d <- subset(d, tpm>0)
		ctss.cluster.df <- .paraclu3(ctss.df = d, minStability = minStability, maxLength = maxLength, removeSingletons = removeSingletons, keepSingletonsAbove = keepSingletonsAbove, reduceToNonoverlapping = reduceToNonoverlapping, useMulticore = useMulticore, nrCores = nrCores)
		ctss.cluster.list[[s]] <- ctss.cluster.df
		
	}
	
	return(ctss.cluster.list)
	
}


#########################################################################
# Imposing predefined clusters to segment the data attached to sequences


.cluster.ctss.strand.predef <- function(ctss.df, custom.clusters) {
	
	clusters <- data.frame(cluster = c(1:length(custom.clusters$start)), chr = custom.clusters$chr, start = custom.clusters$start, end = custom.clusters$end, strand = custom.clusters$strand)
	clusters.ir <- IRanges(start = clusters$start+1, end = clusters$end)
	ctss.df <- ctss.df[order(ctss.df$pos),]
	ctss.df.ir <- IRanges(start = ctss.df$pos, end = ctss.df$pos)
	o <- findOverlaps(clusters.ir, ctss.df.ir)
	ctss.df.1 <- cbind(ctss.df[subjectHits(o),], cluster = clusters$cluster[queryHits(o)], start = clusters$start[queryHits(o)], end = clusters$end[queryHits(o)])
	ctss.df.2 <- data.frame(chr = clusters$chr[!(clusters$cluster %in% clusters$cluster[queryHits(o)])], pos = rep(NA, length(clusters$cluster[!(clusters$cluster %in% clusters$cluster[queryHits(o)])])), strand = clusters$strand[!(clusters$cluster %in% clusters$cluster[queryHits(o)])], tpm = rep(0, length(clusters$cluster[!(clusters$cluster %in% clusters$cluster[queryHits(o)])])), cluster = clusters$cluster[!(clusters$cluster %in% clusters$cluster[queryHits(o)])], start = clusters$start[!(clusters$cluster %in% clusters$cluster[queryHits(o)])], end = clusters$end[!(clusters$cluster %in% clusters$cluster[queryHits(o)])])
	ctss.df <- rbind(ctss.df.1, ctss.df.2)
	ctss.df <- ctss.df[order(ctss.df$cluster),]
	invisible(gc())
	return(ctss.df)
	
} 

.cluster.ctss.chr.predef <- function(ctss.df, custom.clusters) {
	
	ctss.df.p <- subset(ctss.df, strand == "+")
	custom.clusters.p <- subset(custom.clusters, strand == "+")
	if(nrow(custom.clusters.p) > 0) {
		ctss.df.plus <- .cluster.ctss.strand.predef(ctss.df = ctss.df.p, custom.clusters = custom.clusters.p)
		last.plus.cluster <- max(ctss.df.plus$cluster)
	}else{
		ctss.df.plus <- data.frame()
		last.plus.cluster <- 0
	}
	invisible(gc())
	ctss.df.m <- subset(ctss.df, strand == "-")
	custom.clusters.m <- subset(custom.clusters, strand == "-")
	if(nrow(custom.clusters.m) > 0) {			
		ctss.df.minus <- .cluster.ctss.strand.predef(ctss.df = ctss.df.m, custom.clusters = custom.clusters.m)
		ctss.df.minus$cluster <- last.plus.cluster + ctss.df.minus$cluster
	}else{
		ctss.df.minus <- data.frame()
	}
	invisible(gc())
	return(rbind(ctss.df.plus, ctss.df.minus))
	
}

.ctss2clusters.predef <- function(ctss.df, custom.clusters, useMulticore = FALSE, nrCores = NULL) {
	
	if(useMulticore == TRUE){
		library(parallel)
		if(is.null(nrCores)){
			nrCores <- detectCores()
		}		
		ctss.cluster.list <- mclapply(as.list(unique(custom.clusters$chr)), function(x) {
									  
									  ctss.df.chr <- subset(ctss.df, chr == x)
									  custom.clusters <- subset(custom.clusters, chr == x)
									  if(nrow(custom.clusters)>0){
										ctss.cluster.chr.df <- .cluster.ctss.chr.predef(ctss.df = ctss.df.chr, custom.clusters = custom.clusters)
									  }
									  }, mc.cores = nrCores
									  )		
	}else{
		ctss.cluster.list <- lapply(as.list(unique(custom.clusters$chr)), function(x) {
									
									ctss.df.chr <- subset(ctss.df, chr == x)
									custom.clusters <- subset(custom.clusters, chr == x)
									if(nrow(custom.clusters)>0){
										ctss.cluster.chr.df <- .cluster.ctss.chr.predef(ctss.df = ctss.df.chr, custom.clusters = custom.clusters)
									}
									}
									)
	}
	max.clust <- unlist(lapply(ctss.cluster.list, function(x) {max(x$cluster)}))
	max.clust <- cumsum(c(0, max.clust[-length(max.clust)]))
	invisible(gc())
	
	ctss.cluster.list <- lapply(as.list(1:length(max.clust)), function(x) {
								
								a = ctss.cluster.list[[x]]
								a$cluster = a$cluster + max.clust[x]
								return(a)
								
								}
								)
	
	ctss.cluster.df <- do.call(rbind, ctss.cluster.list)
	invisible(gc())
	return(ctss.cluster.df)
	
}	

.summarize.clusters.predef <- function(ctss.cluster.df) {
	
	ctss.cluster <- data.table(ctss.cluster.df)
	ctss.cluster <- ctss.cluster[, list(chr[1], start[1], end[1], strand[1], length(pos), pos[which(tpm == max(tpm))[ceiling(length(which(tpm == max(tpm)))/2)]], sum(tpm), max(tpm)), by = cluster]
	setnames(ctss.cluster, c("cluster", "chr", "start", "end", "strand", "nr_ctss", "dominant_ctss", "tpm", "tpm.dominant_ctss")) 
	ctss.cluster <- data.frame(ctss.cluster)
	ctss.cluster$nr_ctss[ctss.cluster$tpm == 0] <- 0
	
	invisible(gc())
	return(ctss.cluster)
	
}

.predefined.clusters <- function(data, sample.labels, custom.clusters, useMulticore = FALSE, nrCores = NULL){

	ctss.cluster.list <- list()
	for(s in sample.labels) {
		
		message("\t->  ", s)
		d <- data[,c("chr", "pos", "strand", s)]
		colnames(d) <- c("chr", "pos", "strand", "tpm")
		d <- subset(d, tpm > 0)
		ctss.cluster.df <- .ctss2clusters.predef(ctss.df = d, custom.clusters = custom.clusters, useMulticore = useMulticore, nrCores = nrCores)
		ctss.cluster.list[[s]] <- ctss.cluster.df
		invisible(gc())
		
	}
	
	if(useMulticore == TRUE){
		ctss.cluster.list <- mclapply(ctss.cluster.list, .summarize.clusters.predef, mc.cores = nrCores)
	}else{
		ctss.cluster.list <- lapply(ctss.cluster.list, .summarize.clusters.predef)
	}
	invisible(gc())
	return(ctss.cluster.list)
	
}

#' @name tagClusterConvertors
#' 
#' @title Private functions to convert TC formats
#' 
#' @details Interconvert tag clusters (TC) formats used in classes CAGEset
#' (\code{data.frame}) and CAGEexp (\code{GRanges}).
#' 
#' @family df2granges converters
#' 
#' @examples 
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' df <- CAGEr:::getTagCluster(object)[[1]]
#' head(df)
#' gr <- TCdataframe2granges(df)
#' gr
#' head(TCgranges2dataframe(gr))
#' # No round-trip because start and end were not integer in df.
#' if (! identical(df, TCgranges2dataframe(gr))) stop("No round-trip between TCdataframe2granges and TCgranges2dataframe")
#' if (! all(df == TCgranges2dataframe(gr))) stop("No round-trip between TCdataframe2granges and TCgranges2dataframe")
NULL

#' @name CCgranges2dataframe
#' 
#' @rdname tagClusterConvertors
#' 
#' @param gr Consensus clusters in \code{GRanges} format.

TCgranges2dataframe <- function(gr) {
  gr$score <- NULL
  if (!is.null(gr$cluster)) {
    df <- data.frame(cluster = gr$cluster)
    gr$cluster <- NULL
  } else {
    df <- data.frame(cluster = names(gr))
  }
  df <- cbind(df , data.frame( chr     = decode(seqnames(gr))
                             , start   = start(gr)
                             , end     = end(gr)
                             , strand  = decode(droplevels(strand(gr)))))
  df <- cbind(df, mcols(gr))
  as.data.frame(df)
}

#' @name CCdataframe2granges
#' @rdname tagClusterConvertors
#' 
#' @param df Consensus clusters in \code{data.frame} format.

TCdataframe2granges <- function(df) {
	gr <- GRanges( seqnames           = df$chr
	             , ranges             = IRanges(df$start, df$end)
 	             , score              = df$tpm
               , strand             = df$strand)
	mcols(gr) <- cbind( mcols(gr)
	                  , df[,setdiff(colnames(df), c("chr", "start", "end", "strand")), drop = FALSE])
	names(gr) <- rownames(df)
	gr
}
