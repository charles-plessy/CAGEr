#####
# Function that calculates reversed cumulative sums given a list of cumulative sums
# ARGUMENTS: cumsum.list - list of Rle vectors (IRanges package) with cumulative sums (first number in the vector needs to be a zero) (such as returned by 'get.cumsum' function)
# RETURNS: list of Rle vectors containing reversed cumulative sums for all elements in the input cumsum list (Rle vectors are shorter by 1 than original vectors because first zero (i.e. here last number) is omitted) 

.reverse.cumsum <- function(cumsum.list, use.multicore = F, nrCores = NULL) {
	
	if(use.multicore){
		library(multicore)
		if(is.null(nrCores)){
			nrCores <- multicore:::detectCores()
		}		
		cumsum.reversed <- mclapply(cumsum.list, function(x) {cumsum(rev(x[2:length(x)] - x[1:(length(x)-1)]))}, mc.cores = nrCores)
	}else{
		cumsum.reversed <- lapply(cumsum.list, function(x) {cumsum(rev(x[2:length(x)] - x[1:(length(x)-1)]))})
	}
	return(cumsum.reversed)
	
}


.get.dominant.ctss <- function(v, isCumulative = FALSE){
	if(all(unique(v) == 0)){
		idx <- NA
	}else{
		if(isCumulative){
			raw <- v[2:length(v)] - v[1:(length(v)-1)]
		}else{
			raw <- v
		}
	idx <- which(raw == max(raw))
	if(length(idx > 1)){
		idx <- idx[ceiling(length(idx)/2)]
	}
	}
	return(as.integer(idx))
}


.score.promoter.shifting <- function(cum.sum.stages, use.multicore = F, nrCores = NULL) {	
	
	if(use.multicore){
		library(multicore)
		if(is.null(nrCores)){
			nrCores <- multicore:::detectCores()
		}				
		scores <- unlist(mclapply(cum.sum.stages, function(x) {less.tpm <- which(x[nrow(x),] == min(x[nrow(x),]))[1]; if(max(x[,less.tpm])>0) {max(x[,less.tpm] - x[,(3-less.tpm)])/max(x[,less.tpm])}else{return(as.numeric(NA))}}, mc.cores = nrCores))
	}else{
		scores <- unlist(lapply(cum.sum.stages, function(x) {less.tpm <- which(x[nrow(x),] == min(x[nrow(x),]))[1]; if(max(x[,less.tpm])>0) {max(x[,less.tpm] - x[,(3-less.tpm)])/max(x[,less.tpm])}else{return(as.numeric(NA))}}))	
	}

}	
	

