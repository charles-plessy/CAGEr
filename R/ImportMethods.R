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
			
	sample.labels <- sampleLabels(object)
	names(sample.labels) <- rainbow(n = length(sample.labels))

	if(object@inputFilesType == "bam") {
		
		reference.genome <- genomeName(object)
		
		if(reference.genome %in% rownames(installed.packages()) == FALSE){
			stop("Requested genome is not installed! Please install required BSgenome package before running CAGEr.")
		}else if(!paste("package:", reference.genome, sep = "") %in% search()){
			stop("Requested genome is not loaded! Load the genome by calling 'library(", reference.genome, ")'")
		}else{
			genome <- get(ls(paste("package:", reference.genome, sep="")))
		}

		
		bam.files <- inputFiles(object)

		for(f in bam.files) {
			if(file.exists(f) == FALSE){
				stop("Could not locate input file ", f)
			}				
		}
					
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
					CTSS <- .remove.added.G(reads.GRanges.plus, reads.GRanges.minus, genome, correctSystematicG = correctSystematicG)
			}else{
		
				CTSS.plus <- data.frame(chr = as.character(seqnames(reads.GRanges.plus)), pos = as.integer(start(reads.GRanges.plus)), strand = rep("+", times = length(reads.GRanges.plus)), stringsAsFactors = F)
				CTSS.minus <- data.frame(chr = as.character(seqnames(reads.GRanges.minus)), pos = as.integer(end(reads.GRanges.minus)), strand = rep("-", times = length(reads.GRanges.minus)), stringsAsFactors = F)
				CTSS <- rbind(CTSS.plus, CTSS.minus)
				CTSS$tag_count <- 1
				CTSS <- data.table(CTSS)
				CTSS <- CTSS[, as.integer(sum(tag_count)), by = list(chr, pos, strand)]
			
			}
		
			setnames(CTSS, c("chr", "pos", "strand", sample.labels[i])) 
			setkey(CTSS, chr, pos, strand)

			message("\t-> Making CTSSs and counting number of tags...")
				
			library.sizes <- c(library.sizes, as.integer(sum(data.frame(CTSS)[,4])))
		
			if(first == TRUE) {
				CTSS.all.samples <- CTSS
			}else{
				CTSS.all.samples <- merge(CTSS.all.samples, CTSS, all.x = TRUE, all.y = TRUE)
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
				
		for(i in 1:length(ctss.files)) {
			
			message("\nReading in file: ", ctss.files[i], "...")
			
			CTSS <- read.table(file = ctss.files[i], header = F, sep = "\t", colClasses = c("character", "integer", "character", "integer"), col.names = c("chr", "pos", "strand", sample.labels[i]))
			CTSS <- data.table(CTSS)
			setkeyv(CTSS, cols = c("chr", "pos", "strand"))
			
			if(first == TRUE) {
				CTSS.all.samples <- CTSS
				first <- FALSE
			}else{
				CTSS.all.samples <- merge(CTSS.all.samples, CTSS, all.x = TRUE, all.y = TRUE)
			}			
			
		}
		
		CTSS.all.samples <- data.frame(CTSS.all.samples)
		library.sizes <- as.integer(colSums(CTSS.all.samples[,c(4:ncol(CTSS.all.samples)), drop = F], na.rm = T))
		
	}else if(object@inputFilesType == "CTSStable"){
	
		ctss.table.file <- inputFiles(object)
		
		if(length(ctss.table.file) > 1){
			stop("Only one file should be provided when inputFilesType = \"CTSStable\"!")
		}
		if(file.exists(ctss.table.file) == FALSE){
			stop("Could not locate input file ", ctss.table.file)
		}			
		
		CTSS.all.samples <- read.table(file = ctss.table.file, header = F, stringsAsFactors = FALSE)
		if(ncol(CTSS.all.samples) != (length(sample.labels) + 3)){
			stop("Number of provided sample labels must match the number of samples in the CTSS table!")
		}
		library.sizes <- as.integer(apply(CTSS.all.samples[,c(4:ncol(CTSS.all.samples)),drop=F], 2, sum))

	}else{
		
		stop("'inputFilesType' must be one of the supported file types (\"bam\", \"ctss\", \"CTSStable\")")
		
	}
	
	CTSS.all.samples <- data.frame(CTSS.all.samples)
	for(i in 4:ncol(CTSS.all.samples)){
		CTSS.all.samples[is.na(CTSS.all.samples[,i]),i] <- as.integer(0)
	}
	colnames(CTSS.all.samples) <- c("chr", "pos", "strand", sample.labels)
	CTSS.all.samples <- CTSS.all.samples[order(CTSS.all.samples$chr, CTSS.all.samples$pos),]
	rownames(CTSS.all.samples) <- c(1:nrow(CTSS.all.samples))
	
	names(library.sizes) <- sample.labels
	object@sampleLabels <- sample.labels
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
def=function(source, dataset, group, sample){
	standardGeneric("importPublicData")
}
)

setMethod("importPublicData",
signature(source = "character", dataset = "character", sample = "character"),
function (source, dataset, group, sample){
	
	if(source == "ENCODE"){
		
		if("ENCODEprojectCAGE" %in% rownames(installed.packages()) == FALSE){
			stop("Requested CAGE data package is not installed! Please install and load the ENCODEprojectCAGE package, which is available for download from http://promshift.genereg.net/CAGEr/PackageSource/.")
		}else if(!("package:ENCODEprojectCAGE" %in% search())){
			stop("Requested CAGE data package is not loaded! Please load the data package by calling 'library(ENCODEprojectCAGE)'")
		}
			
		if(dataset[1] == "ENCODEtissueCAGEfly"){
			genome.name <- "BSgenome.Dmelanogaster.UCSC.dm3"
			data(ENCODEtissueCAGEfly, envir = environment())
			if(group == "embryo"){
				if(sample == "mixed_embryos_0-24hr"){
					
					ctssTable <- ENCODEtissueCAGEfly[["embryo"]]
					
				}else{
					stop("Specified sample not valid! The dataset 'ENCODEtissueCAGEfly' containes only one sample named 'mixed_embryos_0-24hr'!")
				}
			}else{
				stop("Specified group not valid! The dataset 'ENCODEtissueCAGEfly' containes only one group named 'embryo'!")
			}
		
		}else{
			genome.name <- "BSgenome.Hsapiens.UCSC.hg19"
			data(ENCODEhumanCellLinesSamples, envir = environment())
			info.df <- ENCODEhumanCellLinesSamples
			
			if(!(all(dataset %in% info.df$dataset))){
				stop("Specified dataset(s) not found! Call data(ENCODEhumanCellLinesSamples) and check 'dataset' column for available ENCODE datasets!")
			}
			if(length(dataset) == 1){
				
				if(!(all(group %in% info.df[info.df$dataset == dataset,"group"]))){
					stop("Some of the provided groups cannot be found in the specified dataset!")
				}
				if(length(group) == 1){
					if(!(all(sample %in% info.df[info.df$group == group,"sample"]))){
						stop("Some of the provided samples cannot be found in the specified group!")
					}
				}else if(length(group) == length(sample)){
					if(!all(sapply(c(1:length(group)), function(x) {sample[x] %in% info.df[info.df$group == group[x],"sample"]}))){
						stop("Provided 'group' and 'sample' do not match! Some of the provided samples cannot be found in the corresponding groups!")
					}
				}else{
					stop("Number of elements in the 'group' must be either 1 or must match the number of elements in the 'sample'!")
				}
								
								
			}else if(length(dataset) == length(group)){
				if(!all(sapply(c(1:length(dataset)), function(x) {group[x] %in% info.df[info.df$dataset == dataset[x],"group"]}))){
					stop("Provided 'dataset' and 'group' do not match! Some of the provided groups cannot be found in the corresponding datasets!")
				}
				if(length(group) == length(sample)){
					if(!all(sapply(c(1:length(group)), function(x) {sample[x] %in% info.df[info.df$group == group[x],"sample"]}))){
						stop("Provided 'group' and 'sample' do not match! Some of the provided samples cannot be found in the corresponding groups!")
					}
				}else{
					stop("Number of elements in the 'group' must match the number of elements in the 'sample'!")
				}
			
				
			}else{
				stop("Number of elements in the 'dataset' must be either 1 or must match the number of elements in the 'group' and in the 'sample'!")
			}
			
			data(list = dataset, envir = environment())
			
			if(length(unique(dataset))>1){
				
				for(i in 1:length(unique(dataset))){
					dset <- get(unique(dataset)[i])
					for(j in 1:length(unique(group[which(dataset == unique(dataset)[i])]))){
						g <- unique(group[which(dataset == unique(dataset)[i])])[j]
						ctss <- dset[[g]][, c("chr", "pos", "strand", sample[which((group == g) & (dataset == unique(dataset)[i]))])]
						discard <- apply(ctss[, c(4:ncol(ctss)), drop = F], 1, function(x) {sum(x > 0)}) >= 1
						ctss <- data.table(ctss[discard,])
						setkeyv(ctss, cols = c("chr", "pos", "strand"))
						if(i == 1 & j == 1){
							ctssTable <- ctss							
						}else{
							ctssTable <- merge(ctssTable, ctss, all.x = T, all.y = T)
						}
					}
				}
				ctssTable[is.na(ctssTable)] <- 0
				
			}else{
				
				selected.dataset <- get(dataset)
				if(length(unique(group))>1){
					
					for(i in 1:length(unique(group))){
						g <- unique(group)[i]
						ctss <- selected.dataset[[g]][, c("chr", "pos", "strand", sample[which(group == g)])]
						discard <- apply(ctss[, c(4:ncol(ctss)), drop = F], 1, function(x) {sum(x > 0)}) >= 1
						ctss <- data.table(ctss[discard,])
						setkeyv(ctss, cols = c("chr", "pos", "strand"))
						if(i == 1){
							ctssTable <- ctss							
						}else{
							ctssTable <- merge(ctssTable, ctss, all.x = T, all.y = T)
						}						
					}
					
					ctssTable[is.na(ctssTable)] <- 0					
					
				}else{
					ctssTable <- selected.dataset[[group]][,c("chr", "pos", "strand", sample)]
				}
			}
			
			ctssTable <- data.frame(ctssTable, stringsAsFactors = F, check.names = F)
		}
		
		
	}else if(source == "FANTOM3and4"){
		
		if("FANTOM3and4CAGE" %in% rownames(installed.packages()) == FALSE){
			stop("Requested CAGE data package is not installed! Please install and load the FANTOM3and4CAGE package available from Bioconductor.")
		}else if(!("package:FANTOM3and4CAGE" %in% search())){
			stop("Requested CAGE data package is not loaded! Please load the data package by calling 'library(FANTOM3and4CAGE)'")
		}
		
		data(FANTOMhumanSamples, envir = environment())
		info.df1 <- FANTOMhumanSamples
		data(FANTOMmouseSamples, envir = environment())
		info.df2 <- FANTOMmouseSamples
				
		if(!(all(dataset %in% info.df1$dataset) | all(dataset %in% info.df2$dataset))){
			stop("Specified dataset(s) not found! Call data(FANTOMhumanSamples) and data(FANTOMmouseSamples) and check 'dataset' column for available ENCODE datasets!")
		}
		if(length(grep("human", dataset))>0){
			genome.name <- "BSgenome.Hsapiens.UCSC.hg18"
			info.df <- info.df1
		}else if(length(grep("mouse", dataset))>0){
			genome.name <- "BSgenome.Mmusculus.UCSC.mm9"
			info.df <- info.df2
		}
		
		if(length(dataset) == 1){
			
			if(!(all(group %in% info.df[info.df$dataset == dataset,"group"]))){
				stop("Some of the provided groups cannot be found in the specified dataset!")
			}
			if(length(group) == 1){
				if(!(all(sample %in% info.df[info.df$group == group,"sample"]))){
					stop("Some of the provided samples cannot be found in the specified group!")
				}
			}else if(length(group) == length(sample)){
				if(!all(sapply(c(1:length(group)), function(x) {sample[x] %in% info.df[info.df$group == group[x],"sample"]}))){
					stop("Provided 'group' and 'sample' do not match! Some of the provided samples cannot be found in the corresponding groups!")
				}
			}else{
				stop("Number of elements in the 'group' must be either 1 or must match the number of elements in the 'sample'!")
			}
			
			
		}else if(length(dataset) == length(group)){
			if(!all(sapply(c(1:length(dataset)), function(x) {group[x] %in% info.df[info.df$dataset == dataset[x],"group"]}))){
				stop("Provided 'dataset' and 'group' do not match! Some of the provided groups cannot be found in the corresponding datasets!")
			}
			if(length(group) == length(sample)){
				if(!all(sapply(c(1:length(group)), function(x) {sample[x] %in% info.df[info.df$group == group[x],"sample"]}))){
					stop("Provided 'group' and 'sample' do not match! Some of the provided samples cannot be found in the corresponding groups!")
				}
			}else{
				stop("Number of elements in the 'group' must match the number of elements in the 'sample'!")
			}
			
			
		}else{
			stop("Number of elements in the 'dataset' must be either 1 or must match the number of elements in the 'group' and in the 'sample'!")
		}
		
		data(list = dataset, envir = environment())
		
		if(length(unique(dataset))>1){
			
			for(i in 1:length(unique(dataset))){
				dset <- get(unique(dataset)[i])
				for(j in 1:length(unique(group[which(dataset == unique(dataset)[i])]))){
					g <- unique(group[which(dataset == unique(dataset)[i])])[j]
					ctss <- dset[[g]][, c("chr", "pos", "strand", sample[which(group == g)])]
					discard <- apply(ctss[, c(4:ncol(ctss)), drop = F], 1, function(x) {sum(x > 0)}) >= 1
					ctss <- data.table(ctss[discard,])
					setnames(ctss, c("chr", "pos", "strand", paste(g, sample[which(group == g)], sep = "__")))
					setkeyv(ctss, cols = c("chr", "pos", "strand"))
					if(i == 1 & j == 1){
						ctssTable <- ctss							
					}else{
						ctssTable <- merge(ctssTable, ctss, all.x = T, all.y = T)
					}
				}
			}
			
			ctssTable[is.na(ctssTable)] <- 0
			
		}else{
			
			selected.dataset <- get(dataset)
			if(length(unique(group))>1){
				
				for(i in 1:length(unique(group))){
					g <- unique(group)[i]
					ctss <- selected.dataset[[g]][, c("chr", "pos", "strand", sample[which(group == g)])]
					discard <- apply(ctss[, c(4:ncol(ctss)), drop = F], 1, function(x) {sum(x > 0)}) >= 1
					ctss <- data.table(ctss[discard,])
					setnames(ctss, c("chr", "pos", "strand", paste(g, sample[which(group == g)], sep = "__")))
					setkeyv(ctss, cols = c("chr", "pos", "strand"))
					if(i == 1){
						ctssTable <- ctss							
					}else{
						ctssTable <- merge(ctssTable, ctss, all.x = T, all.y = T)
					}						
				}
				
				ctssTable[is.na(ctssTable)] <- 0					
				
			}else{
				ctssTable <- selected.dataset[[group]][,c("chr", "pos", "strand", sample)]
				setnames(ctssTable, c("chr", "pos", "strand", paste(group, sample, sep = "__")))
			}
		}
		
		ctssTable <- data.frame(ctssTable, stringsAsFactors = F, check.names = F)
	
			
		
	}else if (source == "FANTOM5"){
		
		if(length(dataset) != 1){
			stop("For FANTOM5 only one dataset can be specified and it can be either 'human' or 'mouse'!")
		}else if(!(dataset %in% c("human", "mouse"))){
			stop("For FANTOM5, dataset can be either 'human' or 'mouse'!")
		}
		if(dataset == "human"){
			data(FANTOM5humanSamples, envir = environment())
			samples.info <- FANTOM5humanSamples
			genome.name <- "BSgenome.Hsapiens.UCSC.hg19"
		}else if(dataset == "mouse"){
			data(FANTOM5mouseSamples, envir = environment())
			samples.info <- FANTOM5mouseSamples
			genome.name <- "BSgenome.Mmusculus.UCSC.mm9"
		}		

		if(!(all(sample %in% samples.info$sample))){
			stop(paste("Some sample names cannot be found for the specified dataset! Call data(FANTOM5", dataset, "Samples) and check the 'sample' column for valid sample names!", sep = ""))
		}
		
		
		
		for(i in c(1:length(sample))){
			
			message("Fetching sample: ", sample[i], "...")
			sample.url <- samples.info[samples.info$sample == sample[i], "data_url"]
			con <- gzcon(url(paste(sample.url)))
			ctss <- scan(con, what = list(character(), NULL, integer(), NULL, integer(), character()))
			ctss.df <- data.table(chr = ctss[[1]], pos = ctss[[3]], strand = ctss[[6]], tagCount = ctss[[5]])
			setnames(ctss.df, c("chr", "pos", "strand", sample[i]))
			setkeyv(ctss.df, cols = c("chr", "pos", "strand"))
			if(i == 1){
				ctss.table <- ctss.df
			}else{
				message("Adding sample to CTSS table...\n")
				ctss.table <- merge(ctss.table, ctss.df, all.x = T, all.y = T)
				ctss.table[is.na(ctss.table)] <- 0
			}
			
		}
		
		ctssTable <- data.frame(ctss.table, stringsAsFactors = F, check.names = F)
		
		
	}else if (source == "ZebrafishDevelopment"){

		if("ZebrafishDevelopmentalCAGE" %in% rownames(installed.packages()) == FALSE){
			stop("Requested CAGE data package is not installed! Please install and load the ZebrafishDevelopmentalCAGE package, which is available for download from http://promshift.genereg.net/CAGEr/PackageSource/.")
		}else if(!("package:ZebrafishDevelopmentalCAGE" %in% search())){
			stop("Requested CAGE data package is not loaded! Please load the data package by calling 'library(ZebrafishDevelopmentalCAGE)'")
		}
		
		data(ZebrafishSamples, envir = environment())
		if(dataset == "ZebrafishCAGE"){
			if(group == "development"){
				if(!(all(sample %in% ZebrafishSamples$sample))){
					stop("Some sample names cannot be found for the specified dataset! Call data(ZebrafishSamples) and check the 'sample' column for valid sample names!")
				}else{
					genome.name <- "BSgenome.Drerio.UCSC.danRer7"
					data(ZebrafishCAGE, envir = environment())
					ctssTable <- ZebrafishCAGE[["development"]][,c("chr", "pos", "strand", sample)]
                    ctssTable <- ctssTable[apply(ctssTable[,4:ncol(ctssTable),drop=FALSE], 1, function(x) {any(x>0)}),]
				}
			}else{
				stop("Invalid group name! There is only one group in this dataset named 'development'.")
			}
		}else{
			stop("Invalid dataset name! There is only one available dataset named 'ZebrafishCAGE'.")
		}
		
		
	}else{
		stop("Currently only the following public CAGE data resources are supported: 'FANTOM5', 'FANTOM3and4', 'ENCODE', 'ZebrafishDevelopment'. Refer to CAGEr vignette on how to use those resources!")
	}
	
    rownames(ctssTable) <- c(1:nrow(ctssTable))
    
	sample.labels <- colnames(ctssTable)[4:ncol(ctssTable)]
	names(sample.labels) <- rainbow(n = length(sample.labels))
	myCAGEset <- new("CAGEset", genomeName = genome.name, inputFiles = paste(source, sample.labels, sep = "__"), inputFilesType = source, sampleLabels = sample.labels)
	myCAGEset@librarySizes <- as.integer(colSums(ctssTable[,4:ncol(ctssTable),drop=FALSE]))
	myCAGEset@CTSScoordinates <- ctssTable[, c("chr", "pos", "strand")]
	myCAGEset@tagCountMatrix <- ctssTable[,4:ncol(ctssTable),drop=FALSE]
	
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
	}else if((unique(substr(colors, start = 1, stop = 1)) == "#") & all(unique(unlist(strsplit(substr(colors, start = 2, stop = sapply(colors, width)), split = ""))) %in% c(seq(0,9,1), "A", "B", "C", "D", "E", "F"))){
		names(sample.labels) <- colors
	}else{
		stop("'colors' argument must be a vector of valid color names in R or a vector of hexadecimal specifications (e.g. #008F0AFF). See colors() for a complete list of valid color names.")
	}
	
	object@sampleLabels <- sample.labels
	assign(objName, object, envir = parent.frame())
	invisible(1)
	
}
)


