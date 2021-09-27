#' @include Multicore.R

#' @name scoreShift
#' 
#' @title Calculate promoter shifting score
#' 
#' @description Calculates the shifting score for all consensus clusters (promoters) between
#' two specified (groups of) CAGE datasets.  Shifting score is a measure of differential
#' usage of TSSs within consensus cluster between two samples, which indicates the degree of
#' physical separation of TSSs used in these samples within given consensus cluster. In
#' addition to shifting score, a statistical significance (P-value and FDR) of differential
#' TSS usage is calculated for each consensus cluster using Kolmogorov-Smirnov test.
#' 
#' @param object A \code{\link{CAGEr}} object.
#' 
#' @param groupX,groupY Character vector of the one or more CAGE dataset labels in the first
#' (\code{groupX}) and in the second group (\code{groupY}).  Shifting score for each consensus
#' cluster will be calculated by comparing CAGE signal in the samples from \code{groupX} against
#' the signal in the samples from \code{groupY}.  If there is more than one CAGE dataset in the
#' group, the datasets within that group will be merged together before comparison with the
#' other group.  See Details.
#' 
#' @param testKS Logical, should Kolomogorov-Smirnov test for statistical significance of
#' differential TSS usage be performed, and P-values and FDR returned.  See Details.
#' 
#' @param useTpmKS Logical, should normalized (tpm) values (\code{TRUE}) or raw tag counts
#' (\code{FALSE}) be used to derive sample sizes for Kolomogorov-Smirnov test.  Used only when
#' \code{testKS = TRUE}, otherwise ignored.  See Details.
#' 
#' @param useMulticore Logical, should multicore be used.  \code{useMulticore = TRUE} is
#' supported only on Unix-like platforms.
#' 
#' @param nrCores Number of cores to use when \code{useMulticore = TRUE}.  Default value
#' \code{NULL} uses all detected cores.
#' 
#' @details TSSs within one consensus cluster (promoter) can be used differently in different
#' samples (cell types, tissues, developmental stages), with respect to their position and
#' frequency of usage detected by CAGE.  This function calculates shifting scores of all
#' consensus clusters between two specified (groups of) CAGE samples to detect promoters that
#' are used differently in these two samples.  Shifting score is a measure of differential
#' TSS usage defined as:
#' 
#' \code{score = max(F1 - F2) / max(F1)}
#' 
#' where F1 is a cumulative sum of CAGE signal along consensus cluster in the group of samples
#' with lower total signal in that consensus cluster, and F2 in the opposite group.  Since
#' cumulative sum can be calculated in both forward (5' -> 3') and reverse (3' -> 5')
#' direction, shifting score is calculated for both cases and the bigger value is selected as
#' final shifting score.  Value of the shifting score is in the range \code{[-Inf, 1]}, where
#' value of \code{1} means complete physical separation of TSSs used in the two samples for
#' given consensus cluster.  In general, any non-negative value of the shifting score can be
#' interpreted as the proportion of transcription initiation in the sample with lower expression
#' that is happening "outside" (either upstream or downstream) of the region used for
#' transcription initiation in the other sample.  Negative values indicate no physical
#' separation, \emph{i.e.} the region used for transcription initiation in the sample with
#' lower expression is completely contained within the region used for transcription
#' initiation in the other sample.
#' 
#' In addition to shifting score which indicates only physical separation (upstream or
#' downstream shift of TSSs), a more general assessment of differential TSS usage can be
#' obtained by performing a two-sample Kolmogorov-Smirnov test on cumulative sums of CAGE
#' signal along the consensus cluster.  In that case, cumulative sums in both samples are
#' scaled to range `[0,1]` and are considered to be empirical cumulative distribution functions
#' (ECDF) reflecting sampling of TSS positions during transcription initiation.
#' Kolmogorov-Smirnov test is performed to assess whether the two underlying probability
#' distributions differ.  To obtain P-value (\emph{i.e.} the level at which the
#' null-hypothesis can be rejected), sample sizes that generated the ECDFs are required, in
#' addition to actual K-S statistics calculated from ECDFs.  These are derived either from
#' raw tag counts, \emph{i.e.} exact number of times each TSS in the cluster was sampled
#' during sequencing (when \code{useTpmKS = FALSE}), or from normalized tpm values (when
#' \code{useTpmKS = TRUE}). P-values obtained from K-S tests are further adjusted for
#' multiple testing using Benjamini & Hochberg (BH) method and for each P-value a
#' corresponding false-discovery rate (FDR) is also reported.
#' 
#' Since calculation of shifting scores and Kolmogorov-Smirnov test require cumulative sums
#' along consensus clusters, they have to be calculated beforehand by calling
#' \code{\link{cumulativeCTSSdistribution}} function.
#' 
#' The slots \code{shiftingGroupX}, \code{shiftingGroupY} and
#' \code{consensusClustersShiftingScores} of the provided \code{\link{CAGEexp}} object will
#' be occupied by the information on the groups of CAGE datasets that have been compared and
#' shifting scores of all consensus clusters.  Consensus clusters (promoters) with shifting
#' score and/or FDR above specified threshold can be extracted by calling
#' \code{\link{getShiftingPromoters}} function.
#' 
#' @author Vanja Haberle
#' 
#' @seealso \code{\link{cumulativeCTSSdistribution}}
#' @family CAGEr promoter shift functions
#' 
#' @examples
#'# scoreShift( exampleCAGEexp
#'#           , groupX = c("sample1", "sample2")
#'#           , groupY = "sample3"
#'#           , testKS = TRUE, useTpmKS = FALSE)
#'# head(getShiftingPromoters(exampleCAGEexp))
#'
#' @importFrom stats p.adjust
#' @importFrom utils head
#' @export

setGeneric( "scoreShift"
          , function( object, groupX, groupY
                    , testKS = TRUE, useTpmKS = TRUE, useMulticore = F, nrCores = NULL)
              standardGeneric("scoreShift"))

#' @rdname scoreShift

setMethod( "scoreShift", "CAGEexp"
         , function (object, groupX, groupY, testKS, useTpmKS, useMulticore, nrCores) {
	
	objName <- deparse(substitute(object))
	sample.labels <- sampleLabels(object)
	
	if(!(all(groupX %in% sample.labels) & all(groupY %in% sample.labels))){
		stop("'groupX' and 'groupY' must be vectors of sample labels! Check 'sampleLabels(", objName, ")' for labels of your CAGE datasets!")
	}
	
	message("\nCalculating shifting score...")
	a <- object@CTSScumulativesConsensusClusters
	
	# Problem: originally, CTSScumulativesConsensusClusters were in 0-based
	# coordinates, meaning that the first Rle value was always 0.
	# Then I (Charles) rewrote many functions for the introduction of the
	# CAGEexp class and switched to 1-based coordinates (a bit inadvertently
	# at the beginning).  Most CAGEr functions now expect 1-based coordinates,
	# but it is not yet the case of scoreShift(), and I am running out of time.
	# Hence, the ugly hack that follows.
	
	# Are we 0-based ?
	
	firstIsZero <- function(x) decode(head(x,1) == 0)
	allFirstIsZero <- function(x) all(sapply(x, firstIsZero))
	zeroBased <- function(x) all(sapply(x, allFirstIsZero))
	
	# If not, ...
	
	makeZeroBased <- function(l) {
	  lapply(l, function(s) {
	    lapply(s, function(x) {
	      append(0,x)
	    })
	  })
	}
	if(!zeroBased(a)) a <- makeZeroBased(a)
	
	# -- end of ugly hack ---
	
	a <- a[names(a) %in% c(groupX, groupY)]
	b <- object@consensusClusters
	cumsum.list <- bplapply(a, function(x) {
	  n <- names(x)
	  y <- subset(b, !(b$consensus.cluster %in% as.integer(names(x))))
	  if (nrow(y)>0) {
	    nulls <- lapply(as.list(c(1:nrow(y))), function(t) {
	      Rle(rep(0, y[t, "end"] - y[t, "start"] + 1))
	    })
	    x <- append(x, nulls)
	  }
	  names(x) <- c(n, as.character(y$consensus.cluster))
	  return(x)
	}, BPPARAM = CAGEr_Multicore(useMulticore, nrCores))
	
	cumsum.list.r <- bplapply( cumsum.list, .reverse.cumsum
	                         , useMulticore = useMulticore, nrCores = nrCores
	                         , BPPARAM = CAGEr_Multicore(useMulticore, nrCores))
	
	cumsum.matrices.list.f <- bplapply(as.list(names(cumsum.list[[1]])), function(x) {
	  sapply(names(cumsum.list), function(y) {as.numeric(cumsum.list[[y]][[x]])})
	}, BPPARAM = CAGEr_Multicore(useMulticore, nrCores))
	
	cumsum.matrices.list.r <- bplapply(as.list(names(cumsum.list.r[[1]])), function(x) {
	  m <- matrix(ncol = length(names(cumsum.list.r)), sapply(names(cumsum.list.r), function(y) {
	    as.numeric(cumsum.list.r[[y]][[x]])
	  }))
	  colnames(m) <- names(cumsum.list.r)
	  return(m)
	}, BPPARAM = CAGEr_Multicore(useMulticore, nrCores))
	
  .poolGroups <- function(x) {
    if (!is.matrix(x))
      x <- t(x) # Clusters of length 1 produced vectors instead of matrices...
    cbind( groupX = rowSums(x[, groupX, drop = FALSE])
         , groupY = rowSums(x[, groupY, drop = FALSE]))
  }
	
	cumsum.matrices.groups.f <- bplapply( cumsum.matrices.list.f, .poolGroups
	                                    , BPPARAM = CAGEr_Multicore(useMulticore, nrCores))
	
	cumsum.matrices.groups.r <- bplapply( cumsum.matrices.list.r, .poolGroups
	                                    , BPPARAM = CAGEr_Multicore(useMulticore, nrCores))
	
	names(cumsum.matrices.groups.f) <- names(cumsum.list[[1]])
	names(cumsum.matrices.groups.r) <- names(cumsum.list[[1]])
	
	dominant.ctss.pos <- bplapply(as.list(names(cumsum.matrices.groups.f)), function(x) {
	  sapply(c("groupX", "groupY"), function(y) {.get.dominant.ctss(cumsum.matrices.groups.f[[x]][,y], isCumulative = T)})
	}, BPPARAM = CAGEr_Multicore(useMulticore, nrCores))
	
	dominant.ctss.pos <- data.frame(consensus.cluster = as.integer(names(cumsum.matrices.groups.f)), do.call(rbind, dominant.ctss.pos))
	dominant.ctss.pos <- dominant.ctss.pos[order(dominant.ctss.pos$consensus.cluster),]
	colnames(dominant.ctss.pos) <- c("consensus.cluster", "groupX.pos", "groupY.pos")

	clusters.info <- merge(b[,c(1:5)], dominant.ctss.pos, by.x = "consensus.cluster", by.y = "consensus.cluster")
	clusters.info$groupX.pos <- clusters.info$groupX.pos + clusters.info$start
	clusters.info$groupY.pos <- clusters.info$groupY.pos + clusters.info$start

	n <- names(cumsum.matrices.groups.f)
	cumsum.matrices.groups.f <- lapply(cumsum.matrices.groups.f, function(x) {x[-1,,drop=F]})
	
	scores.f <- .score.promoter.shifting(cumsum.matrices.groups.f, useMulticore = useMulticore, nrCores = nrCores)
	scores.r <- .score.promoter.shifting(cumsum.matrices.groups.r, useMulticore = useMulticore, nrCores = nrCores)
	
	scores <- pmax(scores.f, scores.r)
	names(scores) <- n

	groupX.tpm <- unlist(lapply(cumsum.matrices.groups.f, function(x) {max(x[,"groupX"])}))
	groupY.tpm <- unlist(lapply(cumsum.matrices.groups.f, function(x) {max(x[,"groupY"])}))
	scores.df <- data.frame(consensus.cluster = as.integer(names(scores)), shifting.score = scores, groupX.tpm = groupX.tpm, groupY.tpm = groupY.tpm)
	
	clusters.info <- merge(clusters.info, scores.df, by.x = "consensus.cluster", by.y = "consensus.cluster")
	clusters.info <- clusters.info[,c("consensus.cluster", "shifting.score", "groupX.pos", "groupY.pos", "groupX.tpm", "groupY.tpm")]
	
	
	if(testKS){
		
		message("Testing for significance using Kolmogorov-Smirnov test...\n")
		if(useTpmKS){
		
			n <- (clusters.info$groupX.tpm * clusters.info$groupY.tpm)/(clusters.info$groupX.tpm + clusters.info$groupY.tpm)
			names(n) <- clusters.info$consensus.cluster
			
		}else{
			template.tagcount <- rep(0L, length(consensusClustersGR(object)))
			names(template.tagcount) <- consensusClustersGR(object)$consensus.cluster
			tag.count.list <- list()
		
			for(s in c(groupX, groupY)) {
				ctss          <- CTSStagCountGR(object, s)
				ctss          <- ctss[ctss$filteredCTSSidx]
				ctss.clusters <- consensusClustersGR(object, s)
				tag.count     <- .getTotalTagCount(ctss, ctss.clusters)
				tag.count.new <- template.tagcount
				tag.count.new[ctss.clusters$consensus.cluster] <- as.integer(tag.count)  # for some reason at this step the vector is converted to list when run within this function (does not happen when run normally outside the function)!!!
				tag.count.list[[s]] <- unlist(tag.count.new)
				invisible(gc())
			}
			
			tag.count.m <- do.call(cbind, tag.count.list)
			tag.count.m.new <- cbind(groupX = rowSums(tag.count.m[,groupX,drop=F]), groupY = rowSums(tag.count.m[,groupY,drop=F]))
			n <- (tag.count.m.new[,"groupX"] * tag.count.m.new[,"groupY"])/(tag.count.m.new[,"groupX"] + tag.count.m.new[,"groupY"])
			names(n) <- rownames(tag.count.m.new)
		
		}
		ks.stat <- unlist(bplapply(cumsum.matrices.groups.f, function(x) {ks.s <- .ksStat(x)}, BPPARAM = CAGEr_Multicore(useMulticore, nrCores)))
		p.vals <- .ksPvalue(d = ks.stat, n = n[names(cumsum.matrices.groups.f)])
		fdr <- p.adjust(p.vals, method = "BH")
		p.vals <- data.frame(consensus.cluster = as.integer(names(cumsum.matrices.groups.f)), pvalue.KS = p.vals, fdr.KS = fdr)

		clusters.info <- merge(clusters.info, p.vals, by.x = "consensus.cluster", by.y = "consensus.cluster")

	}
	
	
	object@consensusClustersShiftingScores <- clusters.info
	object@shiftingGroupX <- groupX
	object@shiftingGroupY <- groupY
	object
})

#' Select consensus clusters with shifting score above threshold
#' 
#' Extracts consensus clusters with shifting score and/or FDR (adjusted P-value from
#' Kolmogorov-Smirnov test) above specified threshold. Returns their genomic coordinates,
#' total CAGE signal and the position of dominant TSS in the two compared groups of CAGE
#' samples, along with the value of the shifting score, P-value and FDR.  Scores and
#' P-values/FDR have to be calculated beforehand by calling \code{\link{scoreShift}} function.
#' 
#' @param object A \code{\link{CAGEexp}} object.
#' 
#' @param tpmThreshold Consensus clusters with total CAGE signal \code{>= tpmThreshold}
#'        in each of the compared groups will be returned.
#' 
#' @param scoreThreshold Consensus clusters with shifting score \code{>= scoreThreshold}
#'        will be returned. The default value \code{-Inf} returns all consensus clusters
#'        (for which score could be calculated, \emph{i.e.} the ones that have at least
#'        one tag in each of the comapred samples).
#' 
#' @param fdrThreshold Consensus clusters with adjusted P-value (FDR) from
#'        Kolmogorov-Smirnov test \code{>= fdrThreshold} will be returned. The default
#'        value \code{1} returns all consensus clusters (for which K-S test could be
#'        performed, \emph{i.e.} the ones that have at least one tag in each of the
#'        compared samples).
#' 
#' @return Returns a \code{data.frame} of shifting promoters with genomic coordinates and
#' positions of dominant TSS and CAGE signal in the two compared (groups of) samples, along
#' with shifting score and adjusted P-value (FDR).
#' 
#' @author Vanja Haberle
#' 
#' @family CAGEr promoter shift functions
#' 
#' @examples 
#' #head(getShiftingPromoters( exampleCAGEexp, tpmThreshold = 100
#' #                         , scoreThreshold = 0.4, fdrThreshold = 0.01))
#' 
#' @export

setGeneric( "getShiftingPromoters"
          , function( object, tpmThreshold = 0, scoreThreshold = -Inf, fdrThreshold = 1)
              standardGeneric("getShiftingPromoters"))

#' @rdname getShiftingPromoters

setMethod( "getShiftingPromoters", "CAGEexp"
         , function (object, tpmThreshold, scoreThreshold, fdrThreshold) {

	shifting.scores <- object@consensusClustersShiftingScores
	clusters <- object@consensusClusters
	
	sig.shifting <- shifting.scores[ ( shifting.scores$groupX.tpm >= tpmThreshold &
	                                   shifting.scores$groupY.tpm >= tpmThreshold   ) &
	                                  shifting.scores$shifting.score >= scoreThreshold
	                                , ]

	if("fdr.KS" %in% colnames(shifting.scores)){
		sig.shifting <- sig.shifting[sig.shifting$fdr.KS <= fdrThreshold,]
	}
	
	sig.shifting <- merge(clusters[,c(1:5)], sig.shifting, by.x = "consensus.cluster", by.y = "consensus.cluster", all.x = F, all.y = T)
	
	return(sig.shifting)
})
