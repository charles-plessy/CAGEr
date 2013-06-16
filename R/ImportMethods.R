setGeneric(
name="getCTSS",
def=function(object, sequencingQualityThreshold = 10, mappingQualityThreshold = 20, removeFirstG = TRUE, correctSystematicG = TRUE){
	standardGeneric("getCTSS")
}
)

setMethod("getCTSS",
signature(object = "CAGEset"),
function (object, sequencingQualityThreshold = 10, mappingQualityThreshold = 20, removeFirstG = TRUE, correctSystematicG = TRUE){

	if(!is(object,"CAGEset")){
		stop("Need to initialize the CAGEset object")
	}
	
	
	objName <- deparse(substitute(object))
	
	reference.genome <- genomeName(object)
	
	if(paste("package:", reference.genome, sep = "") %in% search()){
		genome <- get(ls(paste("package:", reference.genome, sep="")))
	}else{
		stop("Required genome is not loaded! Load the genome by calling 'library(", reference.genome, ")'")
	}
		
	if(object@inputFilesType == "bam") {
		
	bam.files <- inputFiles(object)

	for(f in bam.files) {
		
		if(file.exists(f) == FALSE){
			stop("Could not locate input file ", f)
		}			
	
	}
		
	sample.labels <- sampleLabels(object)
	names(sample.labels) <- rainbow(n = length(sample.labels))
			
	library.sizes <- vector()
	first <- TRUE
	
	for(i in 1:length(bam.files)) {
	
		message("\nReading in file: ", bam.files[i], "...")
		
		what <- c("rname", "strand", "pos", "qwidth", "seq", "qual", "mapq")
		param <- ScanBamParam(what = what, flag = scanBamFlag(isUnmappedQuery = FALSE))
		bam <- scanBam(bam.files[i], param = param)
		
		message("\t-> Filtering out low quality reads...")
		
		qual <- bam[[1]]$qual
		
		if(length(unique(width(qual)) != 1)){
			uniq.quals <- unique(width(qual))
			quals.list <- lapply(as.list(uniq.quals), function(x) {idx <- width(qual) == x; q.m <- as(qual[idx], "matrix"); q.avg <- as.integer(rowMeans(q.m)); return(list(idx, q.avg))})
			qa.avg <- unlist(lapply(quals.list, function(x) {return(x[[2]])}))
			idx <- unlist(lapply(quals.list, function(x) {return(which(x[[1]]))}))
			qa.avg <- qa.avg[order(idx)]
		}else{
			qa <- as(qual, "matrix")
			qa.avg <- as.integer(rowMeans(qa))
			
		}
		
		reads.GRanges <- GRanges(seqnames = as.vector(bam[[1]]$rname), IRanges(start = bam[[1]]$pos, width = width(bam[[1]]$seq)), strand = bam[[1]]$strand, qual = qa.avg, mapq = bam[[1]]$mapq, seq = bam[[1]]$seq, read.length = width(bam[[1]]$seq))	
		reads.GRanges <- reads.GRanges[seqnames(reads.GRanges) %in% seqnames(genome)]
		reads.GRanges <- reads.GRanges[!(end(reads.GRanges) > seqlengths(genome)[as.character(seqnames(reads.GRanges))])]
		elementMetadata(reads.GRanges)$mapq[is.na(elementMetadata(reads.GRanges)$mapq)] <- Inf
		reads.GRanges.plus <- reads.GRanges[(as.character(strand(reads.GRanges)) == "+" & elementMetadata(reads.GRanges)$qual >= sequencingQualityThreshold) & (as.character(seqnames(reads.GRanges)) %in% seqnames(genome) & elementMetadata(reads.GRanges)$mapq >= mappingQualityThreshold)]
		reads.GRanges.minus <- reads.GRanges[(as.character(strand(reads.GRanges)) == "-" & elementMetadata(reads.GRanges)$qual >= sequencingQualityThreshold) & (as.character(seqnames(reads.GRanges)) %in% seqnames(genome) & elementMetadata(reads.GRanges)$mapq >= mappingQualityThreshold)]
		
		if(removeFirstG == TRUE){
			
			CTSS <- .remove.added.G(reads.GRanges.plus, reads.GRanges.minus, correctSystematicG = correctSystematicG)
		
		}else{
		
			reads.GRanges <- append(reads.GRanges.plus, reads.GRanges.minus)
			CTSS <- data.frame(chr = as.character(seqnames(reads.GRanges)), pos = as.integer(start(reads.GRanges)), strand = as.character(strand(reads.GRanges)), stringsAsFactors = F)
			CTSS$tag_count <- 1
			CTSS <- data.table(CTSS)
			CTSS <- CTSS[, as.integer(sum(tag_count)), by = list(chr, pos, strand)]
			
		}
		
		setnames(CTSS, c("chr", "pos", "strand", sample.labels[i])) 
		setkey(CTSS, chr, pos, strand)

		message("\t-> Making CTSSs and counting number of tags...")
		
		
		library.sizes <- c(library.sizes, as.integer(sum(data.frame(CTSS)[,4])))
		
		if(first == TRUE) {
			CTSS.all.samples <- as.data.frame(CTSS)
		}else{
			CTSS.all.samples <- merge(CTSS.all.samples, as.data.frame(CTSS), by.x = c(1:3), by.y = c(1:3), all.x = TRUE, all.y = TRUE)
		}
		
		first <- FALSE
		
	}
	
	}else if(object@inputFilesType == "ctss") {
	
		first <- TRUE

		ctss.files <- inputFiles(object)
		
		for(f in ctss.files) {
			
			if(file.exists(f) == FALSE){
				stop("Could not locate input file ", f)
			}			
			
		}
		
		sample.labels = sampleLabels(object)
		
		for(i in 1:length(ctss.files)) {
			
			message("\nReading in file: ", ctss.files[i], "...")
			
			CTSS <- read.table(file = ctss.files[i], header = F, sep = "\t", colClasses = c("character", "integer", "character", "integer"), col.names = c("chr", "pos", "strand", sample.labels[i]))
		
			if(first == TRUE) {
				CTSS.all.samples <- CTSS
				first <- FALSE
			}else{
				CTSS.all.samples <- merge(CTSS.all.samples, CTSS, by.x = c(1:3), by.y = c(1:3), all.x = TRUE, all.y = TRUE)
			}			
			
		}
		
		library.sizes <- as.integer(colSums(CTSS.all.samples[,c(4:ncol(CTSS.all.samples)), drop = F], na.rm = T))
		
	}else{
		
		stop("'inputFilesType' must be one of the supported file types (\"bam\", \"ctss\")")
		
	}
	
	CTSS.all.samples <- data.frame(CTSS.all.samples)
	for(i in 4:ncol(CTSS.all.samples)){
		CTSS.all.samples[is.na(CTSS.all.samples[,i]),i] <- as.integer(0)
	}
	CTSS.all.samples <- CTSS.all.samples[order(CTSS.all.samples$chr, CTSS.all.samples$pos),]
	colnames(CTSS.all.samples) <- c("chr", "pos", "strand", sample.labels)
	
	names(library.sizes) <- sample.labels
	object@librarySizes <- library.sizes
	object@CTSScoordinates <- CTSS.all.samples[,c("chr", "pos", "strand")]
	object@tagCountMatrix <- as.data.frame(CTSS.all.samples[,c(4:ncol(CTSS.all.samples)),drop=F])
	
	cat("\n")
	
	assign(objName, object, envir = parent.frame())
	invisible(1)
	
}
)



setGeneric(
name="importPublicData",
def=function(dataset, group, sample){
	standardGeneric("importPublicData")
}
)

setMethod("importPublicData",
signature(dataset = "list", group = "character", sample = "character"),
function (dataset, group, sample){
	
	if(!(all(group %in% names(dataset)))){
		stop("Specified group(s) not found in the provided dataset!")
	}
	if(length(group) == 1){
		if(!(all(sample %in% colnames(dataset[[group]])))){
			stop("Some of the provided samples cannot be found in the specified group!")
		}
		group <- rep(group, times = length(sample))
	}else if(length(group) == length(sample)){
		if(!all(sapply(c(1:length(group)), function(x) {sample[x] %in% colnames(dataset[[group[x]]])}))){
			stop("Provided 'group' and 'sample' do not match! Some of the provided samples cannot be found in the specified groups!")
		}
	}else{
		stop("Number of elements in the 'group' must be either 1 or must match the number of elements in the 'sample'!")
	}
	
	if(length(unique(group)) > 1){
		for(i in 1:length(unique(group))){
			g <- unique(group)[i]
			ctss <- dataset[[g]][, c("chr", "pos", "strand", sample[which(group == g)])]
			discard <- apply(ctss[, c(4:ncol(ctss)), drop = F], 1, function(x) {sum(x > 0)}) >= 1
			ctss <- ctss[discard,]
			colnames(ctss)[4:ncol(ctss)] <- paste(g, colnames(ctss)[4:ncol(ctss)], sep = "__")
			if(i == 1){
				ctssTable <- ctss
			}else{
				ctssTable <- merge(ctssTable, ctss, by.x = c("chr", "pos", "strand"), by.y = c("chr", "pos", "strand"), all.x = T, all.y = T)
			}
			
		}
		ctssTable[is.na(ctssTable)] <- 0
	}else{
		ctssTable <- dataset[[unique(group)]][, c("chr", "pos", "strand", sample)]
		discard <- all(ctssTable[,c(4:ncol(ctss))] > 0)
		ctssTable <- ctssTable[discard,]
	}
	
	datasetName <- deparse(substitute(dataset))
	organismName <- strsplit(datasetName, split = "CAGE")[[1]][2]
	projectName <- substr(datasetName, start = 1, stop = 6)
	genomes <- c("human" = "BSgenome.Hsapiens.UCSC.hg18", "mouse" = "BSgenome.Mmusculus.UCSC.mm9", "fly" = "BSgenome.Dmelanogaster.UCSC.dm3")
	myCAGEset <- new("CAGEset", genomeName = genomes[organismName], inputFiles = paste(projectName, colnames(ctssTable)[4:ncol(ctssTable)], sep = "__"), inputFilesType = projectName, sampleLabels = colnames(ctssTable)[4:ncol(ctssTable)])

	myCAGEset@librarySizes <- as.integer(colSums(ctssTable[,4:ncol(ctssTable)]))
	myCAGEset@CTSScoordinates <- ctssTable[, c("chr", "pos", "strand")]
	myCAGEset@tagCountMatrix <- ctssTable[,4:ncol(ctssTable)]
	
	return(myCAGEset)
	
}
)


setGeneric(
name="setColors",
def=function(object, colors = NULL){
	standardGeneric("setColors")
}
)

setMethod("setColors",
signature(object = "CAGEset"),
function (object, colors = NULL){

	objName <- deparse(substitute(object))
	sample.labels <- sampleLabels(object)
	
	if(length(colors) == 0){
		names(sample.labels) <- rainbow(n = length(sample.labels))
	}else if(length(colors) != length(sample.labels)){
		stop(paste("Number of provided colors must match the number of samples in the CAGEset object, i.e. must be ", length(sample.labels), "!", sep = ""))
	}else if(all(colors %in% colors())){
		rgb.col <- col2rgb(colors)
		names(sample.labels) <- apply(rgb.col, 2, function(x) {rgb(red = x[1], green = x[2], blue = x[3], alpha = 255, maxColorValue = 255)})
	}else if((unique(substr(colors, start = 1, stop = 1)) == "#") & all(unique(unlist(strsplit(substr(colors, start = 2, stop = sapply(cols, width)), split = ""))) %in% c(seq(0,9,1), "A", "B", "C", "D", "E", "F"))){
		names(sample.labels) <- colors
	}else{
		stop("'colors' argument must be a vector of valid color names in R or a vector of hexadecimal specifications (e.g. #008F0AFF). See colors() for a complete list of valid color names.")
	}
	
	object@sampleLabels <- sample.labels
	assign(objName, object, envir = parent.frame())
	invisible(1)
	
}
)


