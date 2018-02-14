#' @name CAGEset-class
#' @docType class
#' @noRd
#' @export
#'
# @import data.table
# @import Rsamtools
# @import S4Vectors
# @import GenomicRanges
# @import IRanges
# Import whole packages above, following earlier versions of CAGEr.  It may be worth
# trying to replace them by ad-hoc importFrom directives.  This would also remove the
# character(0) warning messages when running Roxygen2.
# https://stackoverflow.com/questions/26104401/character0-warnings-when-running-devtoolsload-all-in-rstudio

setClass(Class = "CAGEset",

	representation(
		genomeName = "character",
		inputFiles = "character",
		inputFilesType = "character",
		sampleLabels = "character",
		librarySizes = "integer",
		CTSScoordinates = "data.frame",
		tagCountMatrix = "data.frame",
		normalizedTpmMatrix = "data.frame",
		CTSSexpressionClusteringMethod = "character",
		CTSSexpressionClasses = "character",
		clusteringMethod = "character",
		filteredCTSSidx = "logical",
		tagClusters = "list",
		CTSScumulativesTagClusters = "list",
		tagClustersQuantileLow = "list",
		tagClustersQuantileUp = "list",
		consensusClusters = "data.frame",
		consensusClustersTpmMatrix = "matrix",
		consensusClustersExpressionClusteringMethod = "character",
		consensusClustersExpressionClasses = "character",
		CTSScumulativesConsensusClusters = "list",
		consensusClustersQuantileLow = "list",
		consensusClustersQuantileUp = "list",
		shiftingGroupX = "character",
		shiftingGroupY = "character",
		consensusClustersShiftingScores = "data.frame"		
	),

	prototype(
		genomeName = character(),
		inputFiles = character(),
		inputFilesType = character(),
		sampleLabels = character(),
		librarySizes = integer(),
		CTSScoordinates = data.frame(),
		tagCountMatrix = data.frame(),
		normalizedTpmMatrix = data.frame(),
		CTSSexpressionClusteringMethod = character(),
		CTSSexpressionClasses = character(),
		clusteringMethod = character(),
		filteredCTSSidx = logical(),
		tagClusters = list(),
		CTSScumulativesTagClusters = list(),
		tagClustersQuantileLow = list(),
		tagClustersQuantileUp = list(),
		consensusClusters = data.frame(),
		consensusClustersTpmMatrix = matrix(NA, nrow = 0, ncol = 0),
		consensusClustersExpressionClusteringMethod = character(),
		consensusClustersExpressionClasses = character(),
		CTSScumulativesConsensusClusters = list(),
		consensusClustersQuantileLow = list(),
		consensusClustersQuantileUp = list(),
		shiftingGroupX = character(),
		shiftingGroupY = character(),
		consensusClustersShiftingScores = data.frame()
	),

	validity = function(object) {
#		if(!(object@genomeName %in% suppressWarnings(suppressMessages(BSgenome::installed.genomes()))))
#		return("'genomeName' must be a name of one of the genome packages available in BSgenome! See 'BSgenome::installed.genomes()'")
#		if(object@genomeName %in% rownames(installed.packages()) == FALSE)
#		return("Requested genome is not installed! Please install required BSgenome package before running CAGEr.")
	  supportedTypes <- c("bam", "bamPairedEnd", "bed", "bedmolecule", "ctss", "CTSStable", "FANTOM5", "ENCODE", "FANTOM3and4", "ZebrafishDevelopment")
		if(!(object@inputFilesType %in% supportedTypes))
		return(paste(sQuote("inputFilesType"), "must be one of supported input file types:",
		             paste(sQuote(supportedTypes), collapse = ", "), "."))
		if(!(object@inputFilesType == "CTSStable") & (length(object@sampleLabels) != length(object@inputFiles))) {
			return("Number of provided sample labels must match the number of input files unless inputFilesType = \"CTSStable\"!")
		}
		if(!(all(nzchar(object@sampleLabels))) | !(all(substr(object@sampleLabels, start = 1, stop = 1) %in% c(letters, LETTERS)))){
			stop("All sample labels must be a non-empty strings beginning with a letter!")
		}
		if(length(unique(object@sampleLabels)) != length(object@sampleLabels)) {
			stop("Duplicated sample labels are not allowed!")
		}
		
	}

)

###############################################################
# Function for displaying CAGEset object in user friendly way

#' @name show
#' @noRd
#' @exportMethod show

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
		cat(sapply(sampleLabels(object), function(x) {paste("\t-> ", x, ": ", nrow(tagClusters(object, samples = x)), "\n", sep = "")}))
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

############
# Coercion

#' @name coerce
#' @noRd
#' @exportMethod coerce

setAs("data.frame", "CAGEset",
    function(from){
        
        if(!(ncol(from) >= 3) | !(all(colnames(from)[1:3] == c("chr", "pos", "strand")))){
            stop("First three columns of the input data.frame must contain chromosome name, genomic position and strand of individual TSSs, and must be named 'chr', 'pos' and 'strand', respectively!")
        }
        if(!(ncol(from))>=4){
            stop("Input data.frame needs to contain at least one column with CAGE tag counts, in addition to first three columns specifying chromosome name, genomic position and strand of individual TSSs!")
        }
        if(!(all(nzchar(colnames(from)))) | !(all(substr(colnames(from),
        start = 1, stop = 1) %in% c(letters, LETTERS)))){
            stop("Names of the columns specifying CAGE tag counts in the input data.frame must be non-empty strings beginning with a letter, as they will be used as sample labels in the resulting CAGEset!")
        }
        if(!(is.integer(from[,"pos"]))){
            stop("The 'pos' column in the input data.frame can contain only non-zero integers as these are interpreted as 1-based genomic coordinates of TSSs! Make sure the 'pos' column is of class 'integer'!")
        }else if(any(from[,"pos"] <= 0)){
            stop("The 'pos' column in the input data.frame can contain only non-zero integers as these are interpreted as 1-based genomic coordinates of TSSs!")
        }
        if(!(all(from[,"strand"] %in% c("+", "-")))){
            stop("The 'strand' column in the input data.frame can contain only '+' or '-'!")
        }
        if(!(all(apply(from[,4:ncol(from),drop=FALSE], 2, is.integer)))){
            stop("The columns specifying CAGE tag counts must be non-negative integers! Make sure these columns are of class 'integer'!")
        }else if(any(apply(from[,4:ncol(from),drop=FALSE], 2, function(x) {any(x < 0)}))){
            stop("The columns specifying CAGE tag counts must be non-negative integers!")
        }
        
        ctss.coordinates <- from[, c("chr", "pos", "strand")]
        ctss.coordinates$chr <- as.character(ctss.coordinates$chr)
        ctss.coordinates$pos <- as.integer(ctss.coordinates$pos)
        ctss.coordinates$strand <- as.character(ctss.coordinates$strand)
        sample.labels <- colnames(from)[4:ncol(from)]
        
        myCAGEset <- new("CAGEset", genomeName = "", inputFiles = "data.frame", inputFilesType = "CTSStable", sampleLabels = sample.labels)
        myCAGEset@librarySizes <- as.integer(colSums(from[,4:ncol(from),drop=FALSE]))
        myCAGEset@CTSScoordinates <- ctss.coordinates
        myCAGEset@tagCountMatrix <- from[,4:ncol(from),drop=FALSE]

        return(myCAGEset)

    }
)


######################
# Merge two CAGEsets

#' Merge two CAGEr objects into one
#' 
#' Merges two [`CAGEr`] objects into one by combining the CTSS genomic
#' coordinates and raw tag counts.  The resulting object will contain a union
#' of TSS positions present in the two input objects and raw tag counts for
#' those TSSs in all samples from both input objects.
#' 
#' @param cs1 A `CAGEr` object
#' @param cs2 A `CAGEr` object
#' 
#' @return Note that merging discards all other information present in the
#' two `CAGEr` objects, that is, the merged object will not contain any
#' normalised tag counts, CTSS clusters, quantile positions, etc., so these
#' have to be calculated again by calling the appropriate functions on the
#' merged object.  Also, it is only possible to merge two objects that contain
#' TSS information for the same reference genome and do not share any sample
#' names.
#' 
#' @return Returns a `CAGEset` or `CAGEexp` object, which contains a union of
#' TSS positions present in the two input objects and raw tag counts for those
#' TSSs in all samples from both input objects.
#' 
#' @author Vanja Haberle
#' 
#' @seealso [`CAGEset`], [`CAGEexp`]
#' 
#' @examples 
#' library(BSgenome.Drerio.UCSC.danRer7)
#' 
#' pathsToInputFiles <- system.file("extdata", c("Zf.unfertilized.egg.chr17.ctss",
#'   "Zf.30p.dome.chr17.ctss", "Zf.prim6.rep1.chr17.ctss"), package="CAGEr")
#'   
#' myCAGEset1 <- new("CAGEset", genomeName = "BSgenome.Drerio.UCSC.danRer7",
#' inputFiles = pathsToInputFiles[1:2], inputFilesType = "ctss", sampleLabels =
#' c("sample1", "sample2"))
#' getCTSS(myCAGEset1)
#' 
#' myCAGEset2 <- new("CAGEset", genomeName = "BSgenome.Drerio.UCSC.danRer7",
#' inputFiles = pathsToInputFiles[3], inputFilesType = "ctss", sampleLabels =
#' "sample3")
#' 
#' getCTSS(myCAGEset2)
#' 
#' myCAGEset <- mergeCAGEsets(myCAGEset1, myCAGEset2)
#' 
#' 
#' ce1 <- CAGEexp(genomeName = "BSgenome.Drerio.UCSC.danRer7",
#' inputFiles = pathsToInputFiles[1:2], inputFilesType = "ctss", sampleLabels =
#' c("sample1", "sample2"))
#' getCTSS(ce1)
#' 
#' ce2 <- CAGEexp(genomeName = "BSgenome.Drerio.UCSC.danRer7",
#' inputFiles = pathsToInputFiles[3], inputFilesType = "ctss", sampleLabels =
#' "sample3")
#' 
#' getCTSS(ce2)
#' 
#' ce <- mergeCAGEsets(ce1, ce2)
#' 
#' @export

setGeneric( "mergeCAGEsets"
          , function(cs1, cs2) {
            
  if (genomeName(cs1) != genomeName(cs2))
    stop("Cannot merge two CAGEsets with data from different genomes!")
  
  if(any(sampleLabels(cs1) %in% sampleLabels(cs2)))
    stop("Cannot merge two CAGEsets that share same sample labels!")
  
  standardGeneric("mergeCAGEsets")
})

#' @rdname mergeCAGEsets

setMethod("mergeCAGEsets", 
signature(cs1 = "CAGEset", cs2 = "CAGEset"),
function (cs1, cs2){
    
    ctss1 <- CTSStagCount(cs1)
    ctss2 <- CTSStagCount(cs2)
    
    ctss <- merge(ctss1, ctss2, by.x = c("chr", "pos", "strand"), by.y = c("chr", "pos", "strand"), all.x = T, all.y = T)
    ctss[is.na(ctss)] <- 0
    ctss <- ctss[order(ctss$chr, ctss$pos),]
    for(i in 4:ncol(ctss)){
        ctss[,i] <- as.integer(ctss[,i])
    }
    
    sample.labels <- c(cs1@sampleLabels, cs2@sampleLabels)
    names(sample.labels) <- rainbow(n = length(sample.labels))
    library.sizes <- as.integer(colSums(ctss[,4:ncol(ctss),drop=FALSE]))
    names(library.sizes) <- sample.labels
    myCAGEset <- new("CAGEset", genomeName = cs1@genomeName, inputFiles = "merged CAGEset", inputFilesType = "CTSStable", sampleLabels = sample.labels)
    myCAGEset@librarySizes <- library.sizes
    myCAGEset@CTSScoordinates <- ctss[,c("chr", "pos", "strand")]
    myCAGEset@tagCountMatrix <- ctss[,4:ncol(ctss),drop=FALSE]

    myCAGEset
})

#' @rdname mergeCAGEsets

setMethod( "mergeCAGEsets"
         , signature(cs1 = "CAGEexp", cs2 = "CAGEexp")
         , function (cs1, cs2) {
           
  # First, make a new CTSS SummarizedExperiment.
    
  getListsOfCTSS <- function(object) {
    lapply(sampleList(object), function(name) {
      ctss <- CTSScoordinatesGR(object)
      mcols(ctss) <- NULL
      score(ctss) <- CTSStagCountDF(object)[[name]]
      ctss[score(ctss) != 0]
    })}
  
  l <- GRangesList(c(getListsOfCTSS(cs1), getListsOfCTSS(cs2)))
  
  # Code duplicated from getCTSS
    
  rowRanges <- sort(unique(unlist(l)))
  mcols(rowRanges) <- NULL

  assay <- DataFrame(V1 = Rle(rep(0L, length(rowRanges))))
  
  expandRange <- function(global, local) {
    x <- Rle(rep(0L, length(global)))
    x[global %in% local] <- score(local)
    x
  }
  
  for (i in seq_along(l))
    assay[,i] <- expandRange(rowRanges, l[[i]])
  
  colnames(assay) <- names(l)
  
  se <- SummarizedExperiment( rowRanges = rowRanges
                            , assays    = SimpleList(counts = assay))

  # Second, merge column metadata
  
  df1 <- colData(cs1)
  df2 <- colData(cs2)
  
  commonCols <- intersect(colnames(df1), colnames(df2))
  df <- rbind(df1[,commonCols], df2[,commonCols])
  
  # Then construct a new CAGEexp object
  
  ce <- CAGEexp(colData = df)
  CTSStagCountSE(ce) <- se
  ce
})
