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
		tagClustersInConsensusClusters = "data.frame",
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
		tagClustersInConsensusClusters = data.frame(),
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
#		if(!(object@genomeName %in% suppressWarnings(suppressMessages(BSgenome::available.genomes()))))
#		return("'genomeName' must be a name of one of the genome packages available in BSgenome! See 'available.genomes()'")
#		if(object@genomeName %in% rownames(installed.packages()) == FALSE)
#		return("Requested genome is not installed! Please install required BSgenome package before running CAGEr.")
		if(!(object@inputFilesType %in% c("bam", "ctss", "CTSStable", "FANTOM5", "ENCODE", "FANTOM3and4", "ZebrafishDevelopment")))
		return("'inputFilesType' must be one of supported input file types (\"bam\", \"ctss\", \"CTSStable\")!")
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


############
# Coercion

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


setGeneric(
name="mergeCAGEsets",
def=function(cs1, cs2){
    standardGeneric("mergeCAGEsets")
}
)

setMethod("mergeCAGEsets",
signature(cs1 = "CAGEset", cs2 = "CAGEset"),
function (cs1, cs2){

    if(cs1@genomeName != cs2@genomeName){
        stop("Cannot merge two CAGEsets with data from different genomes!")
    }
    if(any(cs1@sampleLabels %in% cs2@sampleLabels)){
        stop("Cannot merge two CAGEsets that share same sample labels!")
    }
    
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

    return(myCAGEset)
    
}
)




