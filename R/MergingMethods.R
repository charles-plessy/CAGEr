setGeneric(
name="mergeSamples",
def=function(object, mergeIndex, mergedSampleLabels){
	standardGeneric("mergeSamples")
}
)

setMethod("mergeSamples",
signature(object = "CAGEset", mergeIndex = "numeric"),
function (object, mergeIndex, mergedSampleLabels){

	objName <- deparse(substitute(object))
	sample.labels <- sampleLabels(object)
	tag.count <- object@tagCountMatrix
	lib.sizes <- object@librarySizes
	
	if(!(length(mergeIndex) == length(sample.labels))) {
		stop("length of 'mergeIndex' must match number of samples! See 'sampleLabels(\"", objName, "\")' to list your CAGE samples.")
	}
	if(!(length(unique(mergeIndex)) == length(mergedSampleLabels))) {
		stop("numer of provided 'mergedSampleLabels' must match number of unique values provided in 'mergeIndex'!")
	}
	if(min(nchar(sample.labels)) == 0 | !(all(substr(sample.labels, start = 1, stop = 1) %in% c(letters, LETTERS)))){
		stop("'mergedSampleLabels' must contain non-empty strings beginning with a letter!")
	}
	if(length(unique(sample.labels)) != length(sample.labels)) {
		stop("Duplicated sample labels are not allowed!")
	}
	
	mergeIndex <- as.integer(mergeIndex)
	tag.count.matrix.new <- sapply(sort(unique(mergeIndex)), function(x) {cols <- which(mergeIndex == x); a <- rowSums(tag.count[,cols,drop=F]); return(a)})
	lib.sizes.new <- sapply(sort(unique(mergeIndex)), function(x) {cols <- which(mergeIndex == x); a <- sum(lib.sizes[cols]); return(a)})
	names(lib.sizes.new) <- mergedSampleLabels
	colnames(tag.count.matrix.new) <- mergedSampleLabels
	names(mergedSampleLabels) <- rainbow(n = length(mergedSampleLabels))
	
	new.CAGE.set <- suppressWarnings(suppressMessages(new("CAGEset", genomeName = object@genomeName, inputFiles = paste(mergedSampleLabels, "_merged", sep = ""), inputFilesType = object@inputFilesType, sampleLabels = mergedSampleLabels, librarySizes = lib.sizes.new, CTSScoordinates = object@CTSScoordinates, tagCountMatrix = as.data.frame(tag.count.matrix.new))))
	
	assign(objName, new.CAGE.set, envir = parent.frame())
	invisible(1)	
	
}
)

