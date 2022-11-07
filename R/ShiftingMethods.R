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
#' @author Sarvesh Nikumbh
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
	a <- CTSScumulativesCC(object)
	a <- a[names(a) %in% c(groupX, groupY)]
	b <- consensusClustersGR(object)
	
	cumsum.list <- bplapply(a, function(x) {
	  n <- names(x)
	  # y <- subset(b, !(b$consensus.cluster %in% as.integer(names(x))))
	  y <- subset(b, !(names(b) %in% names(x)))
	  
	  if (length(y)>0) {
	    nulls <- lapply(as.list(c(1:length(y))), function(t) {
	      # Rle(rep(0, y[t, "end"] - y[t, "start"] + 1))
	      Rle(rep(0, end(y)[t] - start(y[t] + 1) ))
	    })
	    x <- append(x, nulls)
	  }
	  # names(x) <- c(n, as.character(y$consensus.cluster))
	  # names(x) <- y$consensus.cluster
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
	  sapply(c("groupX", "groupY"), function(y) {
	    .get.dominant.ctss(cumsum.matrices.groups.f[[x]][,y], isCumulative = TRUE)})
	}, BPPARAM = CAGEr_Multicore(useMulticore, nrCores))
	
	# dominant.ctss.pos <- data.frame(consensus.cluster = names(cumsum.matrices.groups.f), do.call(rbind, dominant.ctss.pos))
	dominant.ctss.pos <- data.frame(consensus.cluster = as.data.frame(b)$consensus.cluster, do.call(rbind, dominant.ctss.pos))
	
	## Is this ordering really needed? 
	## With the new names, this ordering does not work!
	# dominant.ctss.pos <- dominant.ctss.pos[order(dominant.ctss.pos$consensus.cluster),]
	colnames(dominant.ctss.pos) <- c("consensus.cluster", "groupX.pos", "groupY.pos")
	
	
	# clusters.info <- merge(b[,c(1:5)], dominant.ctss.pos, by.x = "consensus.cluster", by.y = "consensus.cluster")
	clusters.info <- merge(as.data.frame(b), dominant.ctss.pos, by.x = "consensus.cluster", by.y = "consensus.cluster", sort=FALSE)
	clusters.info$groupX.pos <- clusters.info$groupX.pos + clusters.info$start
	clusters.info$groupY.pos <- clusters.info$groupY.pos + clusters.info$start

	# n <- names(cumsum.matrices.groups.f)
	
	cumsum.matrices.groups.f <- lapply(cumsum.matrices.groups.f, function(x) {x[-1,,drop=F]})
	
	scores.f <- .score.promoter.shifting(cumsum.matrices.groups.f, 
	  useMulticore = useMulticore, nrCores = nrCores)
	scores.r <- .score.promoter.shifting(cumsum.matrices.groups.r, 
	  useMulticore = useMulticore, nrCores = nrCores)
	
	scores <- pmax(scores.f, scores.r)
	names(scores) <- dominant.ctss.pos$consensus.cluster
	
	
	groupX.tpm <- unlist(lapply(cumsum.matrices.groups.f, function(x) {max(x[,"groupX"])}))
	groupY.tpm <- unlist(lapply(cumsum.matrices.groups.f, function(x) {max(x[,"groupY"])}))
	scores.df <- data.frame(consensus.cluster = names(scores), 
	  shifting.score = scores, groupX.tpm = groupX.tpm, groupY.tpm = groupY.tpm)
	
	clusters.info <- merge(clusters.info, scores.df, sort = FALSE, by.x = "consensus.cluster", by.y = "consensus.cluster")
	clusters.info <- clusters.info[,c("consensus.cluster", "shifting.score", "groupX.pos", "groupY.pos", "groupX.tpm", "groupY.tpm")]
	
	
	if(testKS){
	  
		if(useTpmKS){
  		  
  			n <- (clusters.info$groupX.tpm * clusters.info$groupY.tpm)/(clusters.info$groupX.tpm + clusters.info$groupY.tpm)
  			names(n) <- names(cumsum.matrices.groups.f)
			  
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
  				
  				# tag.count.new[ctss.clusters$consensus.cluster] <- 
  				#   as.integer(tag.count)  
  				tag.count.new <- tag.count
  				
  				# for some reason at this step the vector is converted to list when 
  				# run within this function (does not happen when run normally outside the function)!!!
  				tag.count.list[[s]] <- unlist(tag.count.new)
  				invisible(gc())
  			}
  			
  			tag.count.m <- do.call(cbind, tag.count.list)
  			
  			tag.count.m.new <- cbind(groupX = rowSums(tag.count.m[,groupX,drop=F]), 
  			                        groupY = rowSums(tag.count.m[,groupY,drop=F]))
  			n <- (tag.count.m.new[,"groupX"] * tag.count.m.new[,"groupY"])/(tag.count.m.new[,"groupX"] + tag.count.m.new[,"groupY"])
  			names(n) <- names(cumsum.matrices.groups.f)
		}
		
	  
		ks.stat <- unlist(bplapply(cumsum.matrices.groups.f, function(x) {ks.s <- .ksStat(x)}, BPPARAM = CAGEr_Multicore(useMulticore, nrCores)))
		
		p.vals <- .ksPvalue(d = ks.stat, n = n[names(cumsum.matrices.groups.f)])
		fdr <- p.adjust(p.vals, method = "BH")
		
		p.vals <- data.frame(consensus.cluster = clusters.info$consensus.cluster, pvalue.KS = p.vals, fdr.KS = fdr)
		
		clusters.info <- merge(clusters.info, p.vals, by.x = "consensus.cluster", by.y = "consensus.cluster", 
		  sort = FALSE)
	}
	
	
	
	temp_df <- DataFrame(
	  shifting.score = clusters.info$shifting.score,
	  groupX.pos = clusters.info$groupX.pos,
	  groupY.pos = clusters.info$groupY.pos,
	  groupX.tpm = clusters.info$groupX.tpm,
	  groupY.tpm = clusters.info$groupY.tpm)
	
	## adjust all colnames in temp_df at this point
	colnames(temp_df) <- c(paste("shifting.score", groupX, groupY, sep="."),
	  paste(c("groupX", "groupY"), c(groupX, groupY), rep("pos", 2), sep="."),
	  paste(c("groupX", "groupY"), c(groupX, groupY), rep("tpm", 2), sep=".")
	)
	
	if(testKS){
	    temp_df2 <- DataFrame(pvalue.KS = clusters.info$pvalue.KS,
	                          fdr.KS  = clusters.info$fdr.KS)
	    ## adjust colnames for columns in temp_df2
	    colnames(temp_df2) <- paste(c("pvalue.KS", "fdr.KS"), groupX, groupY, sep=".")
	    temp_df <- cbind.DataFrame(temp_df, temp_df2)
	}
  
	prior_rowdata <- rowData(consensusClustersSE(object))
	use_df <- cbind.DataFrame(prior_rowdata, temp_df)
	
	rownames(use_df) <- rownames(rowData(consensusClustersSE(object)))
	rowData(consensusClustersSE(object)) <- use_df
	message("Done")
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
#' @param groupX,groupY Character vector of the one or more CAGE dataset labels in the first
#' (\code{groupX}) and in the second group (\code{groupY}). Shifting promoters for the specified 
#' group pair are returned. 
#' 
#' @param tpmThreshold Consensus clusters with total CAGE signal \code{>= tpmThreshold}
#'        in each of the compared groups will be returned.
#' 
#' @param scoreThreshold Consensus clusters with shifting score \code{>= scoreThreshold}
#'        will be returned. The default value \code{-Inf} returns all consensus clusters
#'        (for which score could be calculated, \emph{i.e.} the ones that have at least
#'        one tag in each of the compared samples).
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
#' @author Sarvesh Nikumbh
#' 
#' @family CAGEr promoter shift functions
#' 
#' @examples 
#' #head(getShiftingPromoters( exampleCAGEexp, tpmThreshold = 100
#' #                         , scoreThreshold = 0.4, fdrThreshold = 0.01))
#' 
#' @export

setGeneric( "getShiftingPromoters"
          , function( object, groupX, groupY, tpmThreshold = 0, scoreThreshold = -Inf, fdrThreshold = 1)
              standardGeneric("getShiftingPromoters"))

#' @rdname getShiftingPromoters

setMethod( "getShiftingPromoters", "CAGEexp"
         , function (object, groupX, groupY, tpmThreshold, scoreThreshold,
           fdrThreshold) {

  
           
  shiftSc_cname <- paste("shifting.score", groupX, groupY, sep=".")
	shifting.scores <- mcols(consensusClustersGR(object))
	
	## check if the group pairs supplied have shifting scores calculated
	if(!shiftSc_cname %in% colnames(shifting.scores)){
	    stop("Shifting score for the supplied groups is not calculated yet. ", 
	      "Please use ", sQuote("scoreShift()"), "first.")
	}
	
	clusters <- rowData(consensusClustersSE(object)) #consensusClustersGR(object)
	gXtpm_cname <- paste("groupX", groupX, "tpm", sep=".")
	gYtpm_cname <- paste("groupY", groupY, "tpm", sep=".")
	
	## Useful to only keep relevant columns in the final data.frame
	sel_cnames <- c("consensus.cluster", "score", "tpm",
	  paste("shifting.score", groupX, groupY, sep="."),
	  paste(c("groupX", "groupY"), c(groupX, groupY), rep("pos", 2), sep="."),
	  paste(c("groupX", "groupY"), c(groupX, groupY), rep("tpm", 2), sep=".")
	)
	
	fdr_cname <- paste("fdr.KS", groupX, groupY, sep=".")
	if(fdr_cname %in% colnames(shifting.scores)){
	  message("Note: P-values and FDR columns available") 
	  sel_cnames <- c(sel_cnames, paste(c("pvalue.KS", "fdr.KS"), groupX, groupY, sep="."))
	}else{
	  message("Note: P-values and FDR columns not available")  
	}
	
	## find which columns are relevant?
	sig.shifting <- 
	  shifting.scores[ ( shifting.scores[, gXtpm_cname] >= tpmThreshold &
	                     shifting.scores[, gYtpm_cname] >= tpmThreshold &
	                     !is.na(shifting.scores[, shiftSc_cname]))
	                                , sel_cnames]

	## leave-out where score is NA
	sig.shifting <- sig.shifting[sig.shifting[, shiftSc_cname] >= scoreThreshold, ]
	
	## 

	if(fdr_cname %in% colnames(shifting.scores)){
	  sig.shifting <- sig.shifting[
	                      !is.na(sig.shifting[, fdr_cname]), sel_cnames] 
		sig.shifting <- sig.shifting[sig.shifting[, fdr_cname] <= fdrThreshold,]
		
	}
	
	## This merging may no longer be needed
	# sig.shifting <- merge(shifting.scores, sig.shifting, by = "consensus.cluster", 
	#   all.x = F, all.y = T)
	# 
	return(sig.shifting)
})
