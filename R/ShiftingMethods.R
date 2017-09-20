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
#' scaled to range [0,1] and are considered to be empirical cumulative distribution functions
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
#' \code{consensusClustersShiftingScores} of the provided \code{\link{CAGEset}} object will
#' be occupied by the information on the groups of CAGE datasets that have been compared and
#' shifting scores of all consensus clusters.  Consensus clusters (promoters) with shifting
#' score and/or FDR above specified threshold can be extracted by calling
#' \code{\link{getShiftingPromoters}} function.
#' 
#' @author Vanja Haberle
#' 
#' @seealso \code{\link{cumulativeCTSSdistribution}}, \code{\link{getShiftingPromoters}}
#' 
#' @examples
#' scoreShift( exampleCAGEset
#'           , groupX = c("sample1", "sample2")
#'           , groupY = "sample3"
#'           , testKS = TRUE, useTpmKS = FALSE)
#' head(getShiftingPromoters(exampleCAGEset))
#'
#' @importFrom parallel detectCores
#' @importFrom stats p.adjust
#' @export

setGeneric( "scoreShift"
          , function( object, groupX, groupY
                    , testKS = TRUE, useTpmKS = TRUE, useMulticore = F, nrCores = NULL)
              standardGeneric("scoreShift"))

setMethod( "scoreShift"
         , signature(object = "CAGEset", groupX = "character", groupY = "character")
         , function (object, groupX, groupY, testKS, useTpmKS, useMulticore, nrCores) {
	
  useMulticore <- .checkMulticore(useMulticore)
	
	objName <- deparse(substitute(object))
	sample.labels <- sampleLabels(object)
	
	if(!(all(groupX %in% sample.labels) & all(groupY %in% sample.labels))){
		stop("'groupX' and 'groupY' must be vectors of sample labels! Check 'sampleLabels(", objName, ")' for labels of your CAGE datasets!")
	}
	
	message("\nCalculating shifting score...")
	a <- object@CTSScumulativesConsensusClusters
	a <- a[names(a) %in% c(groupX, groupY)]
	b <- object@consensusClusters
	if(useMulticore){
		if(is.null(nrCores)){
			nrCores <- detectCores()
		}		
		cumsum.list <- mclapply(a, function(x) {n <- names(x); y <- subset(b, !(consensus.cluster %in% as.integer(names(x)))); if(nrow(y)>0) {nulls <- lapply(as.list(c(1:nrow(y))), function(t) {Rle(rep(0, y[t, "end"] - y[t, "start"] + 1))}); x <- append(x, nulls)}; names(x) <- c(n, as.character(y$consensus.cluster)); return(x)}, mc.cores = nrCores)
	}else{
		cumsum.list <- lapply(a, function(x) {n <- names(x); y <- subset(b, !(consensus.cluster %in% as.integer(names(x)))); if(nrow(y)>0) {nulls <- lapply(as.list(c(1:nrow(y))), function(t) {Rle(rep(0, y[t, "end"] - y[t, "start"] + 1))}); x <- append(x, nulls)}; names(x) <- c(n, as.character(y$consensus.cluster)); return(x)})
	}

	cumsum.list.r <- list()
	for(i in 1:length(cumsum.list)){
		cumsum.list.r[[i]] <- .reverse.cumsum(cumsum.list[[i]], useMulticore = useMulticore, nrCores = nrCores)
	}
	names(cumsum.list.r) <- names(cumsum.list)
	
	if(useMulticore){
		cumsum.matrices.list.f <- mclapply(as.list(names(cumsum.list[[1]])), function(x) {sapply(names(cumsum.list), function(y) {as.numeric(cumsum.list[[y]][[x]])})}, mc.cores = nrCores)
		cumsum.matrices.list.r <- mclapply(as.list(names(cumsum.list.r[[1]])), function(x) {m <- matrix(sapply(names(cumsum.list.r), function(y) {as.numeric(cumsum.list.r[[y]][[x]])}), ncol = length(names(cumsum.list.r))); colnames(m) <- names(cumsum.list.r); return(m)}, mc.cores = nrCores)
		cumsum.matrices.groups.f <- mclapply(cumsum.matrices.list.f, function(x) {cbind(groupX = rowSums(x[,groupX,drop=F]), groupY = rowSums(x[,groupY,drop=F]))}, mc.cores = nrCores)
		cumsum.matrices.groups.r <- mclapply(cumsum.matrices.list.r, function(x) {cbind(groupX = rowSums(x[,groupX,drop=F]), groupY = rowSums(x[,groupY,drop=F]))}, mc.cores = nrCores)
	}else{
		cumsum.matrices.list.f <- lapply(as.list(names(cumsum.list[[1]])), function(x) {sapply(names(cumsum.list), function(y) {as.numeric(cumsum.list[[y]][[x]])})})		
		cumsum.matrices.list.r <- lapply(as.list(names(cumsum.list.r[[1]])), function(x) {m <- matrix(sapply(names(cumsum.list.r), function(y) {as.numeric(cumsum.list.r[[y]][[x]])}), ncol = length(names(cumsum.list.r))); colnames(m) <- names(cumsum.list.r); return(m)})
		cumsum.matrices.groups.f <- lapply(cumsum.matrices.list.f, function(x) {cbind(groupX = rowSums(x[,groupX,drop=F]), groupY = rowSums(x[,groupY,drop=F]))})
		cumsum.matrices.groups.r <- lapply(cumsum.matrices.list.r, function(x) {cbind(groupX = rowSums(x[,groupX,drop=F]), groupY = rowSums(x[,groupY,drop=F]))})
	}
	
	names(cumsum.matrices.groups.f) <- names(cumsum.list[[1]])
	names(cumsum.matrices.groups.r) <- names(cumsum.list[[1]])
	
	if(useMulticore){
		dominant.ctss.pos <- mclapply(as.list(names(cumsum.matrices.groups.f)), function(x) {sapply(c("groupX", "groupY"), function(y) {.get.dominant.ctss(cumsum.matrices.groups.f[[x]][,y], isCumulative = T)})}, mc.cores = nrCores)
	}else{
		dominant.ctss.pos <- lapply(as.list(names(cumsum.matrices.groups.f)), function(x) {sapply(c("groupX", "groupY"), function(y) {.get.dominant.ctss(cumsum.matrices.groups.f[[x]][,y], isCumulative = T)})})	
	}
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
			
			idx <- object@filteredCTSSidx
			ctss.df <- cbind(CTSScoordinates(object)[idx,], object@tagCountMatrix[idx,,drop=F])
		
			ctss.clusters.orig <- merge(object@consensusClusters, object@consensusClustersTpmMatrix, by.x = 1, by.y = 0)
		
			template.tagcount <- as.integer(rep(0, nrow(ctss.clusters.orig)))
			names(template.tagcount) <- ctss.clusters.orig$consensus.cluster
			tag.count.list <- list()
		
			for(s in c(groupX, groupY)) {
			
				d <- ctss.df[,c("chr", "pos", "strand", s)]
				colnames(d) <- c("chr", "pos", "strand", "tagcount")
				d <- subset(d, tagcount>0)
				ctss.clusters <- ctss.clusters.orig[ctss.clusters.orig[,s]>0,]
				tag.count <- .getTotalTagCount(ctss.df = d, ctss.clusters = ctss.clusters)
				tag.count.new <- template.tagcount
				tag.count.new[names(tag.count)] <- as.integer(tag.count)  # for some reason at this step the vector is converted to list when run within this function (does not happen when run normally outside the function)!!!
				tag.count.list[[s]] <- unlist(tag.count.new)
				invisible(gc())
			
			}
			
			tag.count.m <- do.call(cbind, tag.count.list)
			tag.count.m.new <- cbind(groupX = rowSums(tag.count.m[,groupX,drop=F]), groupY = rowSums(tag.count.m[,groupY,drop=F]))
			n <- (tag.count.m.new[,"groupX"] * tag.count.m.new[,"groupY"])/(tag.count.m.new[,"groupX"] + tag.count.m.new[,"groupY"])
			names(n) <- rownames(tag.count.m.new)
		
		}
		
		if(useMulticore){
			ks.stat <- unlist(mclapply(cumsum.matrices.groups.f, function(x) {ks.s <- .ksStat(x)}, mc.cores = nrCores))
		}else{
			ks.stat <- unlist(lapply(cumsum.matrices.groups.f, function(x) {ks.s <- .ksStat(x)}))
		}
		p.vals <- .ksPvalue(d = ks.stat, n = n[names(cumsum.matrices.groups.f)])
		fdr <- p.adjust(p.vals, method = "BH")
		p.vals <- data.frame(consensus.cluster = as.integer(names(cumsum.matrices.groups.f)), pvalue.KS = p.vals, fdr.KS = fdr)

		clusters.info <- merge(clusters.info, p.vals, by.x = "consensus.cluster", by.y = "consensus.cluster")

	}
	
	
	object@consensusClustersShiftingScores <- clusters.info
	object@shiftingGroupX <- groupX
	object@shiftingGroupY <- groupY

	assign(objName, object, envir = parent.frame())
	invisible(1)

}
)

#' getShiftingPromoters
#' @noRd
#' @export

setGeneric(
name="getShiftingPromoters",
def=function(object, tpmThreshold = 0, scoreThreshold = -Inf, fdrThreshold = 1){
	standardGeneric("getShiftingPromoters")
}
)

setMethod("getShiftingPromoters",
signature(object = "CAGEset"),
function (object, tpmThreshold = 0, scoreThreshold = -Inf, fdrThreshold = 1){

	shifting.scores <- object@consensusClustersShiftingScores
	clusters <- object@consensusClusters
	
	sig.shifting <- subset(shifting.scores, (groupX.tpm >= tpmThreshold & groupY.tpm >= tpmThreshold) & shifting.score >= scoreThreshold)

	if("fdr.KS" %in% colnames(shifting.scores)){
		sig.shifting <- subset(sig.shifting, fdr.KS <= fdrThreshold)
	}
	
	sig.shifting <- merge(clusters[,c(1:5)], sig.shifting, by.x = "consensus.cluster", by.y = "consensus.cluster", all.x = F, all.y = T)
	
	return(sig.shifting)
	
}
)


