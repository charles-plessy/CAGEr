#' @include Multicore.R

#' @name scoreShift
#' 
#' @title Calculate promoter shifting score
#' 
#' @description Calculates the shifting score for all _consensus clusters_
#' (promoters) between two specified (groups of) CAGE datasets.  Shifting score
#' is a measure of differential usage of TSSs within consensus cluster between
#' two samples, which indicates the degree of physical separation of TSSs used
#' in these samples within given _consensus cluster_. In addition to shifting
#' score, a statistical significance (P-value and FDR) of differential TSS usage
#' is calculated for each _consensus cluster_ using Kolmogorov-Smirnov test.
#' 
#' @param object A [`CAGEexp`] object.
#' 
#' @param groupX,groupY Character vector of the one or more CAGE dataset labels
#' in the first (`groupX`) and the second group (`groupY`).  Shifting score for
#' each _consensus cluster_ will be calculated by comparing CAGE signal in the
#' samples from `groupX` against the signal in the samples from `groupY`.  If
#' there is more than one CAGE dataset in the group, the datasets within that
#' group will be merged together before comparison with the other group.
#' 
#' @param testKS Logical, should Kolomogorov-Smirnov test for statistical
#' significance of differential TSS usage be performed, and P-values and FDR
#' returned.
#' 
#' @param useTpmKS Logical, should normalized (tpm) values (`TRUE`) or raw tag
#' counts (`FALSE`) be used to derive sample sizes for Kolomogorov-Smirnov test.
#' Used only when `testKS` equals `TRUE`.
#' 
#' @param useMulticore Logical, should multicore be used.  `useMulticore = TRUE`
#' is supported only on Unix-like platforms.
#' 
#' @param nrCores Number of cores to use when `useMulticore = TRUE`.  Default
#' value `NULL` uses all detected cores.
#' 
#' @details TSSs within one _consensus cluster_ (promoter) can be used
#' differently in different samples (cell types, tissues, developmental stages),
#' with respect to their position and frequency of usage detected by CAGE.  This
#' function calculates shifting scores of all consensus clusters between two
#' specified (groups of) CAGE samples to detect promoters that are used
#' differently in these two samples.  Shifting score is a measure of
#' differential TSS usage defined as:
#' 
#' `score = max(F1 - F2) / max(F1)`
#' 
#' where F1 is a cumulative sum of CAGE signal along consensus cluster in the
#' group of samples with lower total signal in that consensus cluster, and F2 in
#' the opposite group.  Since cumulative sum can be calculated in both forward
#' (5' -> 3') and reverse (3' -> 5') direction, shifting score is calculated for
#' both cases and the bigger value is selected as final shifting score.  Value
#' of the shifting score is in the range `[-Inf, 1]`, where value of `1` means
#' complete physical separation of TSSs used in the two samples for given
#' consensus cluster.  In general, any non-negative value of the shifting score\
#' can be interpreted as the proportion of transcription initiation in the
#' sample with lower expression that is happening "outside" (either upstream or
#' downstream) of the region used for transcription initiation in the other
#' sample.  Negative values indicate no physical separation, _i.e._ the region
#' used for transcription initiation in the sample with lower expression is
#' completely contained within the region used for transcription
#' initiation in the other sample.
#' 
#' In addition to shifting score which indicates only physical separation
#' (upstream or downstream shift of TSSs), a more general assessment of
#' differential TSS usage can be obtained by performing a two-sample
#' Kolmogorov-Smirnov test on cumulative sums of CAGE signal along the consensus
#' cluster.  In that case, cumulative sums in both samples are scaled to range
#' `[0,1]` and are considered to be empirical cumulative distribution functions
#' (ECDF) reflecting sampling of TSS positions during transcription initiation.
#' Kolmogorov-Smirnov test is performed to assess whether the two underlying
#' probability distributions differ.  To obtain P-value _i.e._ the level at
#' which the null-hypothesis can be rejected), sample sizes that generated the
#' ECDFs are required, in addition to actual K-S statistics calculated from
#' ECDFs.  These are derived either from raw tag counts, _i.e._ exact number of
#' times each TSS in the cluster was sampled during sequencing (when `useTpmKS`
#' is set to `FALSE`), or from normalized TPM values (when `useTpmKS` is set to
#' `TRUE`). P-values obtained from K-S tests are further adjusted for multiple
#' testing using Benjamini & Hochberg (BH) method and for each P-value a
#' corresponding false-discovery rate (FDR) is also reported.
#' 
#' Since calculation of shifting scores and Kolmogorov-Smirnov test require
#' cumulative sums along consensus clusters, they have to be calculated
#' beforehand by calling [`cumulativeCTSSdistribution`] function.
#' 
#' Consensus clusters (promoters) with shifting score and/or FDR above specified
#' threshold can be extracted by calling [`getShiftingPromoters`] function.
#' 
#' @return A [`CAGEexp`] object in which the `consensusClusters` experiment slot
#' was updated to contain new metadata columns giving information on the
#' _shifting score_, and its _significance_, and the _position_ and
#' _expression value_ of the dominant CTSS for the shifting groups X and Y.  The
#' column names are constructed by pasting sample names to `shifting.score`,
#' `fdr.KS`, `pvalue.KS`, `groupX.pos`, `groupY.pos`, `groupX.tmp` and
#' `groupY.tpm`.  
#' 
#' @author Vanja Haberle
#' @author Sarvesh Nikumbh
#' 
#' @seealso [`cumulativeCTSSdistribution`]
#' @family CAGEr promoter shift functions
#' @family CAGEr object modifiers
#' 
#' @examples
#' scoreShift( exampleCAGEexp
#'           , groupX = c("Zf.unfertilized.egg")
#'           , groupY = "Zf.30p.dome"
#'           , testKS = TRUE, useTpmKS = FALSE)
#'
#' @importFrom stats p.adjust
#' @importFrom utils head
#' @importFrom S4Vectors cbind.DataFrame
#' @export

setGeneric( "scoreShift"
          , function( object, groupX, groupY
                    , testKS = TRUE, useTpmKS = TRUE, useMulticore = F, nrCores = NULL)
              standardGeneric("scoreShift"))

#' @rdname scoreShift
#' @aliases scoreShift,CAGEexp-method
#' @aliases scoreShift,matrix

setMethod( "scoreShift", "CAGEexp"
         , function (object, groupX, groupY, testKS, useTpmKS, useMulticore, nrCores) {
	
	objName <- deparse(substitute(object))
	sample.labels <- sampleLabels(object)
	
	if(!(all(groupX %in% sample.labels) & all(groupY %in% sample.labels))){
		stop("'groupX' and 'groupY' must be vectors of sample labels! Check 'sampleLabels(", objName, ")' for labels of your CAGE datasets!")
	}
	
	message("\nCalculating shifting score...")
	cumsum.list <- CTSScumulativesCC(object)[c(groupX, groupY)]
	ccs.gr <- consensusClustersGR(object)
	
	cumsum.list <- bplapply(cumsum.list, function(x) {
	  # This function checks if come consensus clusters of `ccs.gr` are absent in an
	  # element of `cumsum.list`, in which case it adds them with a cumulative expression
	  # of zero.
	  y <- subset(ccs.gr, !(names(ccs.gr) %in% names(x)))
	  
	  if (length(y)>0) {
	    nulls <- lapply(seq_along(y), function(t) {
	      Rle(rep(0, end(y)[t] - start(y[t] + 1) ))
	    })
	    x <- append(x, nulls)
	  }
	  
	  return(x)
	}, BPPARAM = CAGEr_Multicore(useMulticore, nrCores))
	
	# Reverse the cumulative sums.
	cumsum.list.r <- bplapply( cumsum.list, .reverse.cumsum
	                         , useMulticore = useMulticore, nrCores = nrCores
	                         , BPPARAM = CAGEr_Multicore(useMulticore, nrCores))
	
	# Transform the list of clusters grouped by samples into a
	# list of matrices grouped by clusters
	cumsum.matrices.list.f <- bplapply(names(cumsum.list[[1]]), function(x) {
	  sapply(names(cumsum.list), function(y) {as.numeric(cumsum.list[[y]][[x]])})
	}, BPPARAM = CAGEr_Multicore(useMulticore, nrCores))
	
	
	cumsum.matrices.list.r <- bplapply(names(cumsum.list.r[[1]]), function(x) {
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
	
	dominant.ctss.pos <- bplapply(names(cumsum.matrices.groups.f), function(x) {
	  sapply(c("groupX", "groupY"), function(y) {
	    .get.dominant.ctss(cumsum.matrices.groups.f[[x]][,y], isCumulative = TRUE)})
	}, BPPARAM = CAGEr_Multicore(useMulticore, nrCores))
	
	dominant.ctss.pos <- data.frame(names(ccs.gr), do.call(rbind, dominant.ctss.pos))
	
	## Is this ordering really needed? 
	## With the new names, this ordering does not work!
	# dominant.ctss.pos <- dominant.ctss.pos[order(dominant.ctss.pos$consensus.cluster),]
	
	colnames(dominant.ctss.pos) <- c("consensus.cluster", "groupX.pos", "groupY.pos")
	
	stopifnot(nrow(ccs.gr) == nrow(dominant.ctss.pos))
	clusters.info <- cbind(as.data.frame(ccs.gr), dominant.ctss.pos)
	
	clusters.info$groupX.pos <- clusters.info$groupX.pos + clusters.info$start
	clusters.info$groupY.pos <- clusters.info$groupY.pos + clusters.info$start

	cumsum.matrices.groups.f <- lapply(cumsum.matrices.groups.f, tail, -1)
	
	scores.f <- sapply(cumsum.matrices.groups.f, scoreShift)
	scores.r <- sapply(cumsum.matrices.groups.r, scoreShift)
	
	scores <- pmax(scores.f, scores.r)
	
	groupX.tpm <- sapply(cumsum.matrices.groups.f, function(x) max(x[,"groupX"]))
	groupY.tpm <- sapply(cumsum.matrices.groups.f, function(x) max(x[,"groupY"]))
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
	## if(testKS) ENDS HERE
	
	
	temp_df <- DataFrame(
	  shifting.score = clusters.info$shifting.score,
	  groupX.pos = clusters.info$groupX.pos,
	  groupY.pos = clusters.info$groupY.pos,
	  groupX.tpm = clusters.info$groupX.tpm,
	  groupY.tpm = clusters.info$groupY.tpm)
	
	
	## adjust all colnames in temp_df at this point
	group_x_str <- paste(groupX, collapse="-")
	group_y_str <- paste(groupY, collapse="-")
	colnames(temp_df) <- c(
	  paste("shifting.score", group_x_str, group_y_str, sep="."),
	  paste(c("groupX.pos", "groupY.pos"), c(group_x_str, group_y_str), sep="."),
	  paste(c("groupX.tpm", "groupY.tpm"), c(group_x_str, group_y_str), sep=".")
	)
	
	
	if(testKS){
	    temp_df2 <- DataFrame(pvalue.KS = clusters.info$pvalue.KS,
	                          fdr.KS  = clusters.info$fdr.KS)
	    ## adjust colnames for columns in temp_df2
	    colnames(temp_df2) <- paste(c("pvalue.KS", "fdr.KS"), 
	                    group_x_str, group_y_str, sep=".")
	    temp_df <- cbind.DataFrame(temp_df, temp_df2)
	}
	
	
	prior_rowdata <- rowData(consensusClustersSE(object))
	

	use_df <- cbind.DataFrame(prior_rowdata, temp_df)
  
	## manage duplicate columns
	use_df <- use_df[ , !duplicated(colnames(use_df))]
	##
	
	rownames(use_df) <- rownames(rowData(consensusClustersSE(object)))
	rowData(consensusClustersSE(object)) <- use_df
	
	##
	object
})


#' @rdname scoreShift
#' @aliases scoreShift,matrix
#' 
# Input: a matrix of two columns contaning the cumulative distributions for
# groupX and groupY
setMethod( "scoreShift", "matrix",
           function (object, groupX, groupY, testKS, useTpmKS, useMulticore, nrCores) {
  less.tpm <- which(object[nrow(object),] == min(object[nrow(object),]))[1]
  if (max(object[,less.tpm]) > 0) {
    max(object[,less.tpm] - object[,(3-less.tpm)])/max(object[,less.tpm])
  } else {
    NA
  }
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
#' getShiftingPromoters( exampleCAGEexp
#'                     , groupX = "Zf.unfertilized.egg"
#'                     , groupY = "Zf.30p.dome") |> head()
#' 
#' @export

setGeneric( "getShiftingPromoters"
          , function( object, groupX, groupY, tpmThreshold = 0, scoreThreshold = -Inf, fdrThreshold = 1)
              standardGeneric("getShiftingPromoters"))

#' @rdname getShiftingPromoters

setMethod( "getShiftingPromoters", "CAGEexp"
         , function (object, groupX, groupY, tpmThreshold, scoreThreshold,
           fdrThreshold) {

  group_x_str <- paste(groupX, collapse="-")
  group_y_str <- paste(groupY, collapse="-")
           
  shiftSc_cname <- paste("shifting.score", group_x_str, group_y_str, sep=".")
  
  
  
	shifting.scores <- mcols(consensusClustersGR(object))
	shifting.scores$consensus.cluster <- rownames(shifting.scores)
	
	## check if the group pairs supplied have shifting scores calculated
	if(!shiftSc_cname %in% colnames(shifting.scores)){
	    stop("Shifting score for the supplied groups is not calculated yet. ", 
	      "Please use ", sQuote("scoreShift()"), "first.")
	}
	
	group_x_str <- paste(groupX, collapse="-")
	group_y_str <- paste(groupY, collapse="-")
	
	# clusters <- rowData(consensusClustersSE(object)) #consensusClustersGR(object)
	gXtpm_cname <- paste("groupX.tpm", group_x_str, sep=".")
	gYtpm_cname <- paste("groupY.tpm", group_y_str, sep=".")
	

	## Useful to only keep relevant columns in the final data.frame
	sel_cnames <- c("consensus.cluster", "score", "score",
	  paste("shifting.score", group_x_str, group_y_str, sep="."),
	  paste(c("groupX.pos", "groupY.pos"), c(group_x_str, group_y_str), sep="."),
	  paste(c("groupX.tpm", "groupY.tpm"), c(group_x_str, group_y_str), sep=".")
	)
	
	fdr_cname <- paste("fdr.KS", group_x_str, group_y_str, sep=".")
	if(fdr_cname %in% colnames(shifting.scores)){
	  message("Note: P-values and FDR columns available") 
	  sel_cnames <- c(sel_cnames, paste(c("pvalue.KS", "fdr.KS"), group_x_str
	    , group_y_str, sep="."))
	}else{
	  message("Note: P-values and FDR columns not available")  
	}
	
	
	## find which rows are relevant?
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
