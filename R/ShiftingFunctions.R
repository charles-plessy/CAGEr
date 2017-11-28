#' @include Multicore.R

#####
# Function that calculates reversed cumulative sums given a list of cumulative sums
# ARGUMENTS: cumsum.list - list of Rle vectors (IRanges package) with cumulative sums (first number in the vector needs to be a zero) (such as returned by 'get.cumsum' function)
# RETURNS: list of Rle vectors containing reversed cumulative sums for all elements in the input cumsum list (Rle vectors are shorter by 1 than original vectors because first zero (i.e. here last number) is omitted) 

.uncumsum <- function(x) c(x[1], diff(x))
 
.reverse.cumsum <- function(cumsum.list, useMulticore = F, nrCores = NULL)
  bplapply(cumsum.list, function(x) cumsum(rev(.uncumsum(x)))
          , BPPARAM = CAGEr_Multicore(useMulticore, nrCores))

.get.dominant.ctss <- function(v, isCumulative = FALSE){
  if (all(unique(v) == 0)) return(NA)
  if (isCumulative)
    v <- .uncumsum(v)
  idx <- which(v == max(v))
  if (length(idx > 1))
    idx <- idx[ceiling(length(idx)/2)]
  idx
}

.score.promoter.shifting <- function(cum.sum.stages, useMulticore = F, nrCores = NULL) {
  unlist(bplapply(cum.sum.stages, function(x) {
    less.tpm <- which(x[nrow(x),] == min(x[nrow(x),]))[1]
      if (max(x[,less.tpm]) > 0) {
        max(x[,less.tpm] - x[,(3-less.tpm)])/max(x[,less.tpm])
      } else {
        NA
      }
  }, BPPARAM = CAGEr_Multicore(useMulticore, nrCores)))
}


#####
# Function that calculates total tag count in CAGE clusters
# ARGUMENTS: ctss.df - data frame with one row per CTSS containing at least four columns, *chr (chromosome) *pos (genomic position of CTSSs) *strand (genomic strand) *tagcount (raw CAGE tag count)
#            ctss.clusters - data frame with one row per cluster containing at least 6 columns, *cluster (cluster ID) *chr (chromosome) *start (start position of the cluster) *end (end position of the cluster) *strand (strand) *dominant_ctss (position of dominant peak)
# RETURNS: integer vector of total tag count per cluster 


setGeneric( ".getTotalTagCount"
          , function(ctss, ctss.clusters)
          	  standardGeneric(".getTotalTagCount"))

setMethod( ".getTotalTagCount", "CTSS"
         , function(ctss, ctss.clusters) {
	o <- findOverlaps(ctss.clusters, ctss)
	totalCount <- tapply( decode(score(ctss)[subjectHits(o)])
	                    , ctss.clusters$consensus.cluster[queryHits(o)]
	                    , sum)
	ctss.clusters$total <- 0
	mcols(ctss.clusters)[as.numeric(names(totalCount)),"total"] <- totalCount
	ctss.clusters$total
})

.ksStat <- function(cumsum.matrix){
	z <- cumsum.matrix[,1]/max(cumsum.matrix[,1]) - cumsum.matrix[,2]/max(cumsum.matrix[,2])
	max(abs(z))
}


.pKS2 <- function(x, tol){
	
	k_max <- sqrt(2-log(tol))
	
	for(i in c(1:length(x))){
		
		if(x[i] < 1){
		
			z <- as.double(-1 * pi^2/(8*x[i]^2))
			w <- as.double(log(x[i]))
			k <- seq(1, k_max - 10^(-10), 2)
			s <- as.double(sum(exp(k^2 * z - w)))
			x[i] <- as.double(s*sqrt(2*pi))

		}else{
		
			z <- -2 * x[i]^2
			s <- -1
			k <- 1
			old <- 0
			new <- 1
			while(abs(old - new) > tol){
				old <- new
				new <- new + 2 * s * exp(z * k^2)
				s <- -1 * s
				k <- k + 1
			}
			x[i] <- new
		}
	}
	x
}

.pkstwo <- function(x, tol = 1e-06) {
  if (is.numeric(x)) { 
    x <- as.double(x)
  } else stop("argument 'x' must be numeric")
  p <- rep(0, length(x))
  p[is.na(x)] <- NA
  IND <- which(!is.na(x) & (x > 0))
  if (length(IND) > 0) {
    p[IND] <- .pKS2(x = x[IND], tol = as.double(tol))
  }
  p
}

.ksPvalue <- function(d, n)	1 - .pkstwo(sqrt(n) * d)
