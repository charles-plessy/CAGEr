####################################################################################
# Implementation of simple distance-based clustering of data attached to sequences
#


.cluster.ctss.strand <- function(ctss.df, max.dist) {
	
	ctss.df <- subset(ctss.df, tpm > 0)
	v <- rep(0, max(ctss.df$pos) + max.dist)
	v[ctss.df$pos + max.dist] <- ctss.df$tpm
	v <- Rle(v)
	c.starts <- cumsum(v@lengths)[which(v@lengths >= max.dist & v@values == 0)] - max.dist
	c.ends <- c(cumsum(v@lengths)[(which(v@lengths >= max.dist & v@values == 0) - 1)[-1]], length(v)) - max.dist
	clusters <- data.frame(cluster = c(1:length(c.starts)), start = c.starts, end = c.ends)
	clusters.ir <- IRanges(start = clusters$start, end = clusters$end)
	ctss.df <- ctss.df[order(ctss.df$pos),]
	ctss.df.ir <- IRanges(start = ctss.df$pos - 1, end = ctss.df$pos)
	o <- findOverlaps(clusters.ir, ctss.df.ir)
	ctss.df <- cbind(ctss.df, clusters$cluster[queryHits(o)])
	colnames(ctss.df)[ncol(ctss.df)] <- 'cluster'
	return(ctss.df)
	
} 

.cluster.ctss.chr <- function(ctss.df, max.dist) {
	
	ctss.df.p <- subset(ctss.df, strand == "+")
	if(nrow(ctss.df.p) > 0) {
		ctss.df.plus <- .cluster.ctss.strand(ctss.df = ctss.df.p, max.dist = max.dist)
	}else{
		ctss.df.plus <- data.frame()
	}
	ctss.df.m <- subset(ctss.df, strand == "-")
	if(nrow(ctss.df.m) > 0) {			
		ctss.df.minus <- .cluster.ctss.strand(ctss.df = ctss.df.m, max.dist = max.dist)
		ctss.df.minus$cluster <- max(0, suppressWarnings(max(ctss.df.plus$cluster))) + ctss.df.minus$cluster
	}else{
		ctss.df.minus <- data.frame()
	}
	
	return(rbind(ctss.df.plus, ctss.df.minus))
	
}

.ctss2clusters <- function(ctss.df, max.dist = 20, useMulticore = FALSE, nrCores = NULL) {
	
	if(useMulticore == TRUE){
		library(parallel)
		if(is.null(nrCores)){
			nrCores <- detectCores()
		}		
		ctss.cluster.list <- mclapply(as.list(unique(ctss.df$chr)), function(x) {
									
									ctss.df.chr <- subset(ctss.df, chr == x)
									ctss.cluster.chr.df <- .cluster.ctss.chr(ctss.df = ctss.df.chr, max.dist = max.dist)
								
									}, mc.cores = nrCores
								)		
	}else{
		ctss.cluster.list <- lapply(as.list(unique(ctss.df$chr)), function(x) {
								
									ctss.df.chr <- subset(ctss.df, chr == x)
									ctss.cluster.chr.df <- .cluster.ctss.chr(ctss.df = ctss.df.chr, max.dist = max.dist)
								
									}
								)
	}
	max.clust <- unlist(lapply(ctss.cluster.list, function(x) {max(x$cluster)}))
	max.clust <- cumsum(c(0, max.clust[-length(max.clust)]))
	
	ctss.cluster.list <- lapply(as.list(1:length(max.clust)), function(x) {
								
								a = ctss.cluster.list[[x]]
								a$cluster = a$cluster + max.clust[x]
								return(a)
								
								}
							)
	
	ctss.cluster.df <- do.call(rbind, ctss.cluster.list)
	return(ctss.cluster.df)
	
}	

.summarize.clusters <- function(ctss.cluster.df, removeSingletons = FALSE, keepSingletonsAbove = Inf) {

	ctss.cluster <- data.table(ctss.cluster.df)
	ctss.cluster <- ctss.cluster[, list(chr[1], min(pos)-1, max(pos), strand[1], length(pos), pos[which(tpm == max(tpm))[ceiling(length(which(tpm == max(tpm)))/2)]], sum(tpm), max(tpm)), by = cluster]
	setnames(ctss.cluster, c("cluster", "chr", "start", "end", "strand", "nr_ctss", "dominant_ctss", "tpm", "tpm.dominant_ctss"))
		
	ctss.cluster <- data.frame(ctss.cluster)
	if(removeSingletons){
		ctss.cluster <- subset(ctss.cluster, nr_ctss > 1 | tpm >= keepSingletonsAbove)
		ctss.cluster$cluster <- c(1:nrow(ctss.cluster))
		rownames(ctss.cluster) <- c(1:nrow(ctss.cluster))
	}
	return(ctss.cluster)
	
}


.distclu <- function(data, sample.labels, max.dist = 20, removeSingletons = FALSE, keepSingletonsAbove = Inf, useMulticore = FALSE, nrCores = NULL) {
	
	ctss.cluster.list <- list()
	for(s in sample.labels) {
		
		message("\t-> ", s)
		d <- data[,c("chr", "pos", "strand", s)]
		colnames(d) <- c("chr", "pos", "strand", "tpm")
		d <- subset(d, tpm > 0)
		ctss.cluster.df <- .ctss2clusters(ctss.df = d, max.dist = max.dist, useMulticore = useMulticore, nrCores = nrCores)
		ctss.cluster.list[[s]] <- ctss.cluster.df
		invisible(gc())
		
	}
	
	if(useMulticore == TRUE){
		ctss.cluster.list <- mclapply(ctss.cluster.list, function(x) {.summarize.clusters(ctss.cluster.df = x, removeSingletons = removeSingletons, keepSingletonsAbove = keepSingletonsAbove)}, mc.cores = nrCores)
	}else{
		ctss.cluster.list <- lapply(ctss.cluster.list, function(x) {.summarize.clusters(ctss.cluster.df = x, removeSingletons = removeSingletons, keepSingletonsAbove = keepSingletonsAbove)})
	}
	return(ctss.cluster.list)
}



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
		o <- findOverlaps(clusters.gr, ignoreSelf = TRUE, type = "within")
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

