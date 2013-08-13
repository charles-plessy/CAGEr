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
		if(object@genomeName %in% rownames(installed.packages()) == FALSE)
		return("Requested genome is not installed! Please install required BSgenome package before running CAGEr.")
		if(!(object@inputFilesType %in% c("bam", "ctss", "CTSStable", "FANTOM", "ENCODE")))
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
