################################################################
# Functions for retrieving data from CAGEset and CAGEexp objects

#' @name genomeName
#' 
#' @title Extracting genome name from CAGEr objects
#' 
#' @description Extracts the name of a referent genome from a
#' \code{\link{CAGEset}} and \code{\link{CAGEexp}} objects.
#' 
#' @param object A CAGEset or CAGEexp object.
#' 
#' @return Returns a name of a BSgenome package used as a referent genome.
#' 
#' @family CAGEr accessor methods
#' 
#' @author Vanja Haberle
#' 
#' @examples 
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' genomeName(exampleCAGEset)
#' 
#' @docType methods
#' @rdname genomeName
#' @export

setGeneric(
name="genomeName",
def=function(object){
	standardGeneric("genomeName")
})

setMethod("genomeName",
signature(object = "CAGEset"),
function (object){
	object@genomeName
})

setMethod("genomeName",
signature(object = "CAGEexp"),
function (object){
 metadata(object)$genomeName
})

#' @name inputFiles
#' 
#' @title Extracting paths to input files from CAGEr objects
#' 
#' @description Extracts the paths to CAGE data input files from
#' \code{\link{CAGEset}} and \code{\link{CAGEexp}} objects.
#' 
#' @param object A CAGEset or CAGEexp object.
#' 
#' @return Returns a character vector of paths to CAGE data input files.
#' 
#' @family CAGEr accessor methods
#' 
#' @author Vanja Haberle
#' 
#' @examples 
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' inputFiles(exampleCAGEset)
#' 
#' @docType methods
#' @export

setGeneric(
name="inputFiles",
def=function(object){
	standardGeneric("inputFiles")
})

setMethod("inputFiles",
signature(object = "CAGEset"),
function (object){
	object@inputFiles
})

setMethod("inputFiles",
signature(object = "CAGEexp"),
function (object){
  object$inputFiles
})

#' @name inputFilesType
#'
#' @title Extracting type of input files from CAGEr objects
#' 
#' @description Extracts the information on the type of CAGE data input
#' files from \code{\link{CAGEset}} and \code{\link{CAGEexp}} objects.
#' 
#' @param object A CAGEset or CAGEexp object.
#' 
#' @return Returns the label of the file type of CAGE data input files,
#' \emph{e.g.} \code{"bam"} or \code{"ctss"}.  In the case of \code{CAGEexp}
#' objects, the return value is character vector with one member per sample.
#' 
#' @author Vanja Haberle
#' 
#' @examples 
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' inputFilesType(exampleCAGEset)
#' 
#' @family CAGEr accessor methods
#' @docType methods
#' @export

setGeneric(
name="inputFilesType",
def=function(object){
	standardGeneric("inputFilesType")
})

setMethod("inputFilesType",
signature(object = "CAGEset"),
function (object){
	object@inputFilesType
})

setMethod("inputFilesType",
signature(object = "CAGEexp"),
function (object){
  object$inputFilesType
})

#' @name sampleLabels
#' 
#' @title Extracting CAGE datasets labels from CAGEr objects
#' 
#' @description Extracts the labels of CAGE datasets (samples, experiments)
#' from \code{\link{CAGEset}} and \code{\link{CAGEexp}} objects.
#' 
#' @param object A CAGEset or CAGEexp object.
#' 
#' @return Returns a character vector of labels of all CAGE datasets present
#' in the CAGEr object.
#' 
#' @author Vanja Haberle
#' 
#' @examples 
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' sampleLabels(exampleCAGEset)
#' 
#' @family CAGEr accessor methods
#' @docType methods
#' @export

setGeneric(
name="sampleLabels",
def=function(object){
	standardGeneric("sampleLabels")
})

setMethod("sampleLabels",
signature(object = "CAGEset"),
function (object){
	object@sampleLabels
})

setMethod("sampleLabels",
signature(object = "CAGEexp"),
function (object){
  object$sampleLabels
})

#' @name librarySizes
#' 
#' @title Extracting library sizes from CAGEr objects
#' 
#' @description Extracts the library sizes (total number of CAGE tags) for all CAGE datasets
#' from \code{\link{CAGEset}} and \code{\link{CAGEexp}} objects.
#' 
#' @param object A CAGEset or CAGEexp object.
#' 
#' @return Returns an integer vector of total number of CAGE tags (library size) for all CAGE
#' datasets in the CAGEr object.
#' 
#' @examples 
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' librarySizes(exampleCAGEset)
#'
#' @author Vanja Haberle
#' @family CAGEr accessor methods
#' @docType methods
#' @export

setGeneric(
name="librarySizes",
def=function(object){
	standardGeneric("librarySizes")
})

setMethod("librarySizes",
signature(object = "CAGEset"),
function (object){
	object@librarySizes
})

setMethod("librarySizes",
signature(object = "CAGEexp"),
function (object){
  as.integer(object$librarySizes)
})

#' @name CTSScoordinates
#' 
#' @title Extracting genomic coordinates of TSSs from CAGEr object
#' 
#' @description Extracts the genomic coordinates of all detected TSSs
#' from \code{\link{CAGEset}} and \code{\link{CAGEexp}} objects.
#' 
#' @param object A CAGEset or CAGEexp object.
#' 
#' @return Returns a \code{data.frame} with genomic coordinates of all TSSs. \code{pos}
#' column contains 1-based coordinate of the TSS.
#' 
#' @seealso
#' \code{\link{getCTSS}}
#' 
#' @examples
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' CTSS <- CTSScoordinates(exampleCAGEset)
#' head(CTSS)
#'
#' @author Vanja Haberle
#' @family CAGEr accessor methods
#' @docType methods
#' @export

setGeneric(
name="CTSScoordinates",
def=function(object){
	standardGeneric("CTSScoordinates")
})

setMethod("CTSScoordinates",
signature(object = "CAGEset"),
function (object){
	object@CTSScoordinates
})

setMethod("CTSScoordinates",
signature(object = "CAGEexp"),
function (object){
  gr <- rowRanges(experiments(object)$tagCountMatrix)
  data.frame( chr = as.character(seqnames(gr))
            , pos = start(gr)
            , strand = as.character(strand(gr))
            , stringsAsFactors = FALSE)
})

#' CTSScoordinatesGR
#' 
#' @noRd
#' 
#' Same as CTSScoordinates, but as GRanges
#' 
#' Will be more documented if finally exported

setGeneric(
name="CTSScoordinatesGR",
def=function(object){
	standardGeneric("CTSScoordinatesGR")
})

setMethod("CTSScoordinatesGR",
signature(object = "CAGEset"),
function (object){
  ctssCoord <- object@CTSScoordinates
  ctssCoord <- GRanges(ctssCoord$chr, IRanges(ctssCoord$pos, ctssCoord$pos), ctssCoord$strand)
  genome(ctssCoord) <- object@genomeName
  ctssCoord
})

setMethod("CTSScoordinatesGR",
signature(object = "CAGEexp"),
function (object){
  rowRanges(CTSStagCountSE(object))
})

#' @name CTSStagCount
#' 
#' @title Extracting CAGE tag count for TSSs from CAGEr objects
#' 
#' @description Extracts the tag count for all detected TSSs in all CAGE datasets
#'              from \code{\link{CAGEset}} and \code{\link{CAGEexp}} objects.
#' 
#' @param object A CAGEset or CAGEexp object.
#'  
#' @return Returns a \code{data.frame} with number of CAGE tags supporting each TSS
#' (rows) in every CAGE dataset (columns).
#' 
#' @seealso \code{\link{getCTSS}}
#' 
#' @examples
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' tagCount <- CTSStagCount(exampleCAGEset)
#' head(tagCount)
#' 
#' @author Vanja Haberle
#' @family CAGEr accessor methods
#' @docType methods
#' @export

setGeneric(
name="CTSStagCount",
def=function(object){
	standardGeneric("CTSStagCount")
})

setMethod("CTSStagCount",
signature(object = "CAGEset"),
function (object){
	cbind(object@CTSScoordinates, object@tagCountMatrix)
})

setMethod("CTSStagCount",
signature(object = "CAGEexp"),
function (object){
  cbind( CTSScoordinates(object)
       , as.data.frame(lapply(assay(experiments(object)$tagCountMatrix), as.integer)))
})

#' CTSStagCountDf
#' 
#' @noRd

setGeneric(
name="CTSStagCountDf",
def=function(object){
	standardGeneric("CTSStagCountDf")
})

setMethod("CTSStagCountDf",
signature(object = "CAGEset"),
function (object){
	object@tagCountMatrix
})

setMethod("CTSStagCountDf",
signature(object = "CAGEexp"),
function (object){
  as.data.frame(lapply(assay(experiments(object)$tagCountMatrix), as.integer))
})

#' CTSStagCountDF
#' 
#' Same as CTSStagCountDf, but as DataFrame
#' 
#' @noRd

setGeneric(
name="CTSStagCountDF",
def=function(object){
	standardGeneric("CTSStagCountDF")
})

setMethod("CTSStagCountDF",
signature(object = "CAGEset"),
function (object){
	DF <- object@tagCountMatrix
	DF <- lapply(DF, as.integer)
	DF <- lapply(DF, Rle)
	DataFrame(DF)
})

setMethod("CTSStagCountDF",
signature(object = "CAGEexp"),
function (object){
  assay(CTSStagCountSE(object))
})

#' @name CTSStagCountTable
#' 
#' @title Extracting CAGE tag count for TSSs from CAGEr objects
#' 
#' @description Extracts the tag count for all detected TSSs in all CAGE datasets
#' from \code{\link{CAGEset}} and \code{\link{CAGEexp}} objects.
#' 
#' @param object A CAGEset or CAGEexp object.
#'  
#' @return Returns an expression table with the number of CAGE tags supporting each
#' TSS (rows) in every CAGE dataset (columns).  The table is in \code{data.frame}
#' format for \code{CAGEset} objects and in \code{DataFrame} format for \code{CAGEexp}
#' objects.  Use this function when the next consumer can handle both formats.
#' 
#' @seealso \code{\link{getCTSS}}
#' 
#' @examples
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' tagCount <- CTSStagCount(exampleCAGEset)
#' head(tagCount)
#' 
#' @author Vanja Haberle
#' @family CAGEr accessor methods
#' @docType methods
#' @export

setGeneric(
name="CTSStagCountTable",
def=function(object){
	standardGeneric("CTSStagCountTable")
})

setMethod("CTSStagCountTable",
signature(object = "CAGEset"),
  function(object) CTSStagCountDf(object)
)

setMethod("CTSStagCountTable",
signature(object = "CAGEexp"),
  function(object) CTSStagCountDF(object)
)

#' CTSStagCountSE
#' 
#' Same as CTSStagCount, but as SummarizedExperiment
#' 
#' @noRd

setGeneric(
name="CTSStagCountSE",
def=function(object){
	standardGeneric("CTSStagCountSE")
})

setMethod("CTSStagCountSE",
signature(object = "CAGEset"),
function (object){
	  colData <- data.frame(row.names = sampleLabels(object), samplename = sampleLabels(object), samplecolor = names(sampleLabels(object)))
    SummarizedExperiment(assays = list(counts=CTSStagCountDF(object)), rowData = CTSScoordinatesGR(object), colData = colData)
})

setMethod("CTSStagCountSE",
signature(object = "CAGEexp"),
function (object){
  experiments(object)$tagCountMatrix
})

#' @name CTSSnormalizedTpm
#' 
#' @title Extracting normalized CAGE signal for TSSs from CAGEr objects
#' 
#' @description Extracts the normalized CAGE signal for all detected TSSs
#' in all CAGE datasets from \code{\link{CAGEset}} and
#' \code{\link{CAGEexp}} objects.
#' 
#' @param object A CAGEset or CAGEexp object.
#' 
#' @return Returns a \code{data.frame} with normalized CAGE signal supporting
#' each TSS (rows) in every CAGE dataset (columns).
#' 
#' @seealso \code{\link{normalizeTagCount}}
#' 
#' @examples 
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' CAGEsignal <- CTSSnormalizedTpm(exampleCAGEset)
#' head(CAGEsignal)
#' 
#' @author Vanja Haberle
#' @family CAGEr accessor methods
#' @docType methods
#' @export

setGeneric(
name="CTSSnormalizedTpm",
def=function(object){
	standardGeneric("CTSSnormalizedTpm")
})

setMethod("CTSSnormalizedTpm",
signature(object = "CAGEset"),
function (object){
	cbind(object@CTSScoordinates, object@normalizedTpmMatrix)
})

#' CTSSnormalizedTpmDf
#' 
#' @noRd

setGeneric(
name="CTSSnormalizedTpmDf",
def=function(object){
	standardGeneric("CTSSnormalizedTpmDf")
})

setMethod("CTSSnormalizedTpmDf",
signature(object = "CAGEset"),
function (object){
	object@normalizedTpmMatrix
})

setMethod("CTSSnormalizedTpmDf",
signature(object = "CAGEexp"),
function (object){
  as.data.frame(lapply(assay(experiments(object)$normalizedTpmMatrix), as.integer))
})

#' @name CTSSclusteringMethod
#' 
#' @title Extracting CTSS clustering method from CAGEr objects.
#' 
#' @description Extracts the label of the method used for CTSS clustering into tag
#' clusters from \code{\link{CAGEset}} and \code{\link{CAGEexp}} objects.
#' 
#' @param object A CAGEset or CAGEexp object.
#' 
#' @return Returns a label of the method used for CTSS clustering.
#' 
#' @seealso \code{\link{clusterCTSS}}
#' 
#' @examples 
#' load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))
#' CTSSclusteringMethod(exampleCAGEset)
#' 
#' @author Vanja Haberle
#' @family CAGEr accessor methods
#' @docType methods
#' @export

setGeneric(
name="CTSSclusteringMethod",
def=function(object){
	standardGeneric("CTSSclusteringMethod")
})

setMethod("CTSSclusteringMethod",
signature(object = "CAGEset"),
function (object){
	object@clusteringMethod
})

setMethod("CTSSclusteringMethod",
signature(object = "CAGEexp"),
function (object){
	metadata(object)$clusteringMethod
})

# tagClusters

setGeneric(
name="tagClusters",
def=function(object, sample, returnInterquantileWidth = FALSE, qLow = NULL, qUp = NULL){
	standardGeneric("tagClusters")
})

setMethod("tagClusters",
signature(object = "CAGEset"),
function (object, sample, returnInterquantileWidth = FALSE, qLow = NULL, qUp = NULL){
	if(sample %in% object@sampleLabels){
		tc <- object@tagClusters[[sample]]
        if(returnInterquantileWidth & (length(qLow) == 0 | length(qUp) == 0)){
            stop("No quantiles specified! Please specify which quantile positions should be used to calculate width (qLow and qUp arguments)!")
        }else if(returnInterquantileWidth & (length(object@tagClustersQuantileLow)==0 & length(object@tagClustersQuantileUp)==0)){
            stop("Interquantile width cannot be returned because no quantile positions have been calculated yet! Run 'quantilePositions()' first to get the positions of the desired quantiles!")
		}else if(returnInterquantileWidth & (!(paste("q_", qLow, sep = "") %in% colnames(object@tagClustersQuantileLow[[sample]]) & paste("q_", qUp, sep = "") %in% colnames(object@tagClustersQuantileUp[[sample]])))){
			stop("Interquantile width cannot be returned because specified quantile positions have not been calculated! Run 'quantilePositions()' again to get the positions of the desired quantiles!")
		}else if(returnInterquantileWidth){
			tc.w <- merge(object@tagClustersQuantileLow[[sample]], object@tagClustersQuantileUp[[sample]])
			tc.w <- tc.w[,c(1, which(colnames(tc.w) == paste("q_", qLow, sep = "")), which(colnames(tc.w) == paste("q_", qUp, sep = "")))]
			tc.w$interquantile_width <- tc.w[,3] - tc.w[,2] + 1
			tc <- merge(tc, tc.w)
		}else{
			tc <- object@tagClusters[[sample]]
		}
		return(tc)
	}else{
		stop("Provided 'sample' not in the CAGE set! Check sampleLabels()")
	}
})


setGeneric(
name="consensusClusters",
def=function(object, sample=NULL, returnInterquantileWidth = FALSE, qLow = NULL, qUp = NULL){
    standardGeneric("consensusClusters")
})

setMethod("consensusClusters",
signature(object = "CAGEset"),
function (object, sample = NULL, returnInterquantileWidth = FALSE, qLow = NULL, qUp = NULL){
	
    if(length(sample) == 0){
        return(object@consensusClusters)
    }else if(sample %in% object@sampleLabels){
        
        cc.s <- cbind(cluster = as.integer(rownames(consensusClustersTpm(object))), tpm = consensusClustersTpm(object)[,sample])
        
        if(returnInterquantileWidth & (length(qLow) == 0 | length(qUp) == 0)){
            stop("No quantiles specified! Please specify which quantile positions should be used to calculate width (qLow and qUp arguments)!")
        }else if(returnInterquantileWidth & (length(object@consensusClustersQuantileLow)==0 & length(object@consensusClustersQuantileUp)==0)){
            stop("Interquantile width cannot be returned because no quantile positions for consensus clusters have been calculated yet! Run 'quantilePositions()' first to get the positions of the desired quantiles!")
        }else if(returnInterquantileWidth & (!(paste("q_", qLow, sep = "") %in% colnames(object@consensusClustersQuantileLow[[sample]]) & paste("q_", qUp, sep = "") %in% colnames(object@consensusClustersQuantileUp[[sample]])))){
            stop("Interquantile width cannot be returned because specified quantile positions have not been calculated for consensus clusters! Run 'quantilePositions()' again to get the positions of the desired quantiles!")
        }else if(returnInterquantileWidth){
            cc <- object@consensusClusters

            cc.cumsum <- object@CTSScumulativesConsensusClusters[[sample]]
            a <- lapply(cc.cumsum, function(x) {.get.dominant.ctss(as.numeric(x), isCumulative = T)})
            b <- data.frame(consensus.cluster = as.integer(names(a)), dominant_ctss = unlist(a))
            #cc <- merge(b, cc.s, by.x = 1, by.y = 1, all.x = T, all.y = F)
            
            cc <- merge(cc[,-which(colnames(cc) == "tpm")], b, by.x = 1, by.y = 1, all.x = F, all.y = T)
            cc$dominant_ctss <- cc$start + cc$dominant_ctss
            
            cc <- merge(cc, cc.s, by.x = 1, by.y = 1, all.x = T, all.y = F)
            
            ctss <- CTSSnormalizedTpm(object)[,c("chr", "pos", "strand", sample)]
            cc <- merge(cc, ctss, by.x = c("chr", "strand", "dominant_ctss"), by.y = c("chr", "strand", "pos"), all.x = T, all.y = F)
            colnames(cc)[ncol(cc)] <- "tpm.dominant_ctss"
            cc <- cc[,c("consensus.cluster", "chr", "start", "end", "strand", "dominant_ctss", "tpm", "tpm.dominant_ctss")]
            
            cc.w <- merge(object@consensusClustersQuantileLow[[sample]], object@consensusClustersQuantileUp[[sample]])
            cc.w <- cc.w[,c(1, which(colnames(cc.w) == paste("q_", qLow, sep = "")), which(colnames(cc.w) == paste("q_", qUp, sep = "")))]
            cc.w$interquantile_width <- cc.w[,3] - cc.w[,2] + 1
            cc <- merge(cc, cc.w, by.x = "consensus.cluster", by.y = "cluster", all.x = T)
        }else{
            cc <- object@consensusClusters
            cc <- merge(cc[,-which(colnames(cc) == "tpm")], cc.s, by.x = 1, by.y = 1, all.x = T, all.y = F)
            cc <- subset(cc, tpm>0)
        }
        return(cc)
    }else{
        stop("Provided 'sample' not in the CAGE set! Check sampleLabels()")
    }

})


setGeneric(
name="consensusClustersTpm",
def=function(object){
	standardGeneric("consensusClustersTpm")
})

setMethod("consensusClustersTpm",
signature(object = "CAGEset"),
function (object){
	object@consensusClustersTpmMatrix
})


setGeneric(
name="expressionClasses",
def=function(object, what){
	standardGeneric("expressionClasses")
})

setMethod("expressionClasses",
signature(object = "CAGEset"),
function (object, what){
	if(what == "CTSS"){
		classes <- unique(sort(object@CTSSexpressionClasses))
		if(length(classes)>0){
			return(classes)
		}else{
			stop("No expression clustering of CTSSs has been done yet!")
		}
	}else if(what == "consensusClusters") {
		classes <- unique(sort(object@consensusClustersExpressionClasses))
		if(length(classes)>0){
			return(classes)
		}else{
			stop("No expression clustering of consensus clusters has been done yet!")
		}
	}else{
		stop("'what' parameter must be one of the (\"CTSS\", \"consensusClusters\")")
	}
})

#' GeneExpSE
#' 
#' Retreives the SummarizedExperiment containing gene expression levels.
#' 
#' @noRd for the moment

setGeneric(
name="GeneExpSE",
def=function(object){
	standardGeneric("GeneExpSE")
})

setMethod("GeneExpSE",
signature(object = "CAGEset"),
function (object){
	stop("Not implemented for the CAGEset class.")
})

setMethod("GeneExpSE",
signature(object = "CAGEexp"),
function (object){
  experiments(object)$geneExpMatrix
})

###############################################################
# Function for displaying CAGEset object in user friendly way

setMethod("show", 

	signature(object = "CAGEset"), 
	function(object) {
	cat("\nS4 Object of class CAGEset\n\n")
	cat("=======================================\n")
	cat("Input data information\n")
	cat("=======================================\n")
	cat("Reference genome (organism): ", genomeName(object), "\n", sep = "")
	cat("Input file type: ", inputFilesType(object), "\n", sep = "")
	cat("Input file names: ", paste(inputFiles(object), collapse = ", "), "\n", sep = "")
	cat("Sample labels: ", paste(sampleLabels(object), collapse = ", "), "\n", sep = "")
	cat("=======================================\n")
	cat("CTSS information\n")
	cat("=======================================\n")
	if(nrow(CTSScoordinates(object))>0){
	max_out = min(3, nrow(CTSScoordinates(object)))
		cat("CTSS chromosome: ", paste(paste(CTSScoordinates(object)$chr[1:max_out], collapse = ", "), ", ...", sep = ""), "\n", sep = "")
		cat("CTSS position: ", paste(paste(CTSScoordinates(object)$pos[1:max_out], collapse = ", "), ", ...", sep = ""), "\n", sep = "")
		cat("CTSS strand: ", paste(paste(CTSScoordinates(object)$strand[1:max_out], collapse = ", "), ", ...", sep = ""), "\n", sep = "")
	}else{
		cat("CTSS chromosome:\n")	
		cat("CTSS position:\n")	
		cat("CTSS strand:\n")	
	}
	cat("Tag count:\n")
	if(nrow(object@tagCountMatrix) > 0){
		cat(sapply(c(1:length(sampleLabels(object))), function(x) {paste("\t-> ", sampleLabels(object)[x], ": ", paste(paste(object@tagCountMatrix[1:max_out,x], collapse = ", "), ", ...\n", sep = ""), sep = "")}))
	}
	cat("Normalized tpm:\n")
	if(nrow(object@normalizedTpmMatrix)>0){	
		cat(sapply(c(1:length(sampleLabels(object))), function(x) {paste("\t-> ", sampleLabels(object)[x], ": ", paste(paste(formatC(object@normalizedTpmMatrix[1:max_out,x], format = "f", digits = 3), collapse = ", "), ", ...\n", sep = ""), sep = "")}))
	}
	cat("=======================================\n")
	cat("Tag cluster (TC) information\n")
	cat("=======================================\n")
	cat("CTSS clustering method: ", object@clusteringMethod, "\n", sep = "")
	cat("Number of TCs per sample:\n")
	if(length(object@tagClusters) > 0){
		cat(sapply(sampleLabels(object), function(x) {paste("\t-> ", x, ": ", nrow(tagClusters(object, sample = x)), "\n", sep = "")}))
	}
	cat("=======================================\n")
	cat("Consensus cluster information\n")
	cat("=======================================\n")
	if(nrow(consensusClusters(object))>0){
		max_out = min(3, nrow(consensusClusters(object)))
		cat("Number of consensus clusters: ", nrow(consensusClusters(object)), "\n", sep = "")
		cat("Consensus cluster chromosome: ", paste(paste(consensusClusters(object)$chr[1:max_out], collapse = ", "), ", ...", sep = ""), "\n", sep = "")
		cat("Consensus cluster start: ", paste(paste(consensusClusters(object)$start[1:max_out], collapse = ", "), ", ...", sep = ""), "\n", sep = "")
		cat("Consensus cluster end: ", paste(paste(consensusClusters(object)$end[1:max_out], collapse = ", "), ", ...", sep = ""), "\n", sep = "")
		cat("Consensus cluster strand: ", paste(paste(consensusClusters(object)$strand[1:max_out], collapse = ", "), ", ...", sep = ""), "\n", sep = "")
		cat("Normalized tpm:\n")
		cat(sapply(c(1:length(sampleLabels(object))), function(x) {paste("\t-> ", sampleLabels(object)[x], ": ", paste(paste(formatC(object@consensusClustersTpmMatrix[1:max_out,x], format = "f", digits = 3), collapse = ", "), ", ...\n", sep = ""), sep = "")}))
	}else{
		cat("Number of consensus clusters:\n")
		cat("Consensus cluster chromosome:\n")
		cat("Consensus cluster start:\n")
		cat("Consensus cluster end:\n")
		cat("Consensus cluster strand:\n")
		cat("Normalized tpm:\n")		
	}
	cat("=======================================\n")
	cat("Expression profiling\n")
	cat("=======================================\n")
	cat("Expression clustering method: ", object@consensusClustersExpressionClusteringMethod, "\n", sep = "")
	if(length(object@consensusClustersExpressionClasses) > 0){
		cat("Expression clusters for consensus clusters: ", paste(paste(object@consensusClustersExpressionClasses[1:max_out], collapse = ", "), ", ...", sep = ""), "\n", sep = "")
	}else{
		cat("Expression clusters for consensus clusters:\n")
	}
	cat("=======================================\n")
	cat("Promoter shifting\n")
	cat("=======================================\n")
	if(length(object@shiftingGroupX)>0){
		max_out = min(3, nrow(object@consensusClustersShiftingScores))
		cat("GroupX: ", paste(object@shiftingGroupX, collapse = ", "), "\n", sep = "")
		cat("GroupY: ", paste(object@shiftingGroupY, collapse = ", "), "\n", sep = "")
		cat("Shifting scores: ", paste(formatC(object@consensusClustersShiftingScores$shifting.score[1:max_out], format = "f", digits = 3), collapse = ", "), "\n", sep = "")
		cat("KS p-values (FDR adjusted): ", paste(formatC(object@consensusClustersShiftingScores$fdr.KS[1:max_out], format = "e", digits = 2), collapse = ", "), "\n", sep = "")
	}else{
		cat("GroupX:\n")
		cat("GroupY:\n")
		cat("Shifting scores:\n")
		cat("KS p-values (FDR adjusted):\n")
	}	
	cat("\n")
		
})

