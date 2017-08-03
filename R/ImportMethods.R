#' @include AllClasses.R CAGEexp.R

#' @name getCTSS
#' 
#' @title Reading CAGE data from input file(s) and detecting TSSs
#' 
#' @description Reads input CAGE datasets into CAGEr object, constructs CAGE
#' transcriptions start sites (CTSSs) and counts number of CAGE tags supporting every
#' CTSS in each input experiment.  See \code{\link{inputFilesType}} for details on
#' the supported input formats. Preprocessing and quality filtering of input CAGE
#' tags, as well as correction of CAGE-specific 'G' nucleotide addition bias can be
#' also performed before constructing TSSs.
#' 
#' @param object A \code{\link{CAGEset}} or a A \code{\link{CAGEexp}} object.
#' 
#' @param sequencingQualityThreshold Only CAGE tags with average sequencing quality
#'   \code{>= sequencingQualityThreshold} and mapping quality \code{>=
#'   mappingQualityThreshold} are kept. Used only if \code{inputFileType(object)
#'   == "bam"} or \code{inputFileType(object) == "bamPairedEnd"}, \emph{i.e} when
#'   input files are BAM files of aligned sequenced CAGE tags, otherwise ignored.  If
#'   there are no sequencing quality values in the BAM file (\emph{e.g.} HeliScope
#'   single molecule sequencer does not return sequencing qualities) all reads will by
#'   default have this value set to -1.  Since the default value of
#'   \code{sequencingQualityThreshold} is 10, all the reads will consequently be
#'   discarded. To avoid this behaviour and keep all sequenced reads set
#'   \code{sequencingQualityThreshold} to -1 when processing data without sequencing
#'   qualities.  If there is no information on mapping quality in the BAM file
#'   (\emph{e.g.} software used to align CAGE tags to the referent genome does not
#'   provide mapping quality) the \code{mappingQualityThreshold} parameter is ignored.
#'   In case of paired-end sequencing BAM file (i.e. \code{inputFileType(object) ==
#'   "bamPairedEnd"}) only the first mate of the properly paired reads (i.e. the five
#'   prime end read) will be read and subject to specified thresholds.
#'   
#' @param mappingQualityThreshold See sequencingQualityThreshold.
#' 
#' @param removeFirstG Logical, should the first nucleotide of the CAGE tag be removed
#'   in case it is a G and it does not map to the referent genome (\emph{i.e.} it is a
#'   mismatch).  Used only if \code{inputFileType(object) == "bam"} or
#'   \code{inputFileType(object) == "bamPairedEnd"}, \emph{i.e} when input files are
#'   BAM files of aligned sequenced CAGE tags, otherwise ignored.  See Details.
#' 
#' @param correctSystematicG Logical, should the systematic correction of the first G
#'   nucleotide be performed for the positions where there is a G in the CAGE tag and G
#'   in the genome.  This step is performed in addition to removing the first G of the
#'   CAGE tags when it is a mismatch, \emph{i.e.} this option can only be used when
#'   \code{removeFirstG = TRUE}, otherwise it is ignored.  The frequency of adding a G
#'   to CAGE tags is estimated from mismatch cases and used to systematically correct
#'   the G addition for positions with G in the genome.  Used only if
#'   \code{inputFileType(object) == "bam"} or \code{inputFileType(object) ==
#'   "bamPairedEnd"}, \emph{i.e} when input files are BAM files of aligned sequenced
#'   CAGE tags, otherwise ignored.  See Details.
#' 
#' @details In the CAGE experimental protocol an additional G nucleotide is often attached
#' to the 5' end of the tag by the template-free activity of the reverse transcriptase used
#' to prepare cDNA (Harbers and Carninci, Nature Methods 2005).  In cases where there is a
#' G at the 5' end of the CAGE tag that does not map to the corresponding genome sequence,
#' it can confidently be considered spurious and should be removed from the tag to avoid
#' misannotating actual TSS. Thus, setting \code{removeFirstG = TRUE} is highly recommended.
#' 
#' However, when there is a G both at the beginning of the CAGE tag and in the genome, it is
#' not clear whether the original CAGE tag really starts at this position or the G nucleotide
#' was added later in the experimental protocol.  To systematically correct CAGE tags mapping
#' at such positions, a general frequency of adding a G to CAGE tags can be calculated from
#' mismatch cases and applied to estimate the number of CAGE tags that have G added and
#' should actually start at the next nucleotide/position.  The option \code{correctSystematicG}
#' is an implementation of the correction algorithm described in Carninci \emph{et al.},
#' Nature Genetics 2006, Supplementary Information section 3-e.
#' 
#' @return For \code{\link{CAGEset}} objects, the slots \code{librarySizes}, \code{CTSScoordinates}
#' and \code{tagCountMatrix} will be occupied by the information on CTSSs created from input CAGE
#' files.  For \code{\link{CAGEexp}} objects, the \code{tagCountMatrix} experiment will be
#' occupied by a \code{SummarizedExperiment} containing the expression data as a \code{DataFrame}
#' of \code{Rle} integers, and the CTSS coordinates as a \code{GRanges} object.  In both cases
#' the expression data can be retreived with \code{\link{CTSStagCount*}} functions.  In addition,
#' the library sizes are calculated and stored in the object.
#' 
#' @references 
#' 
#' Harbers and Carninci (2005) Tag-based approaches for transcriptome research and genome
#' annotation, \emph{Nature Methods} \bold{2}(7):495-502.
#' 
#' Carninci \emph{et al.} (2006) Genome-wide analysis of mammalian promoter architecture and
#' evolution, \emph{Nature Genetics} \bold{38}(7):626-635.
#' 
#' @author Vanja Haberle
#' 
#' @seealso \code{\link{CTSScoordinates}}, \code{\link{CTSStagCount}},
#'   \code{\link{CTSStagCountTable}}, \code{\link{inputFilesType}},
#'    \code{\link{librarySizes}}.
#' 
#' @family CAGEr object modifiers
#' 
#' @examples
#' library(BSgenome.Drerio.UCSC.danRer7)
#' 
#' pathsToInputFiles <- system.file("extdata", c("Zf.unfertilized.egg.chr17.ctss",
#'   "Zf.30p.dome.chr17.ctss", "Zf.prim6.rep1.chr17.ctss"), package="CAGEr")
#'   
#' labels <- paste("sample", seq(1,3,1), sep = "")
#' 
#' myCAGEset <- new("CAGEset", genomeName = "BSgenome.Drerio.UCSC.danRer7",
#'  inputFiles = pathsToInputFiles, inputFilesType = "ctss", sampleLabels = labels)
#' 
#' getCTSS(myCAGEset)
#' 
#' @docType methods
#' 
#' @importFrom  data.table setnames
#' @importFrom  data.table setkeyv
#' @importFrom  Rsamtools ScanBamParam
#' @importFrom  Rsamtools scanBamFlag
#' @importFrom  Rsamtools scanBam
#' @importFrom  S4Vectors DataFrame
#' @importFrom  S4Vectors Rle
#' @export

setGeneric(
name="getCTSS",
def=function(object, sequencingQualityThreshold = 10, mappingQualityThreshold = 20, removeFirstG = TRUE, correctSystematicG = TRUE){
	standardGeneric("getCTSS")
}
)

checkRefGenomeIsLoaded <- function(reference.genome) {
  if(reference.genome %in% rownames(installed.packages()) == FALSE){
			stop("Requested genome is not installed! Please install required BSgenome package before running CAGEr.")
		}else if(!paste("package:", reference.genome, sep = "") %in% search()){
			stop("Requested genome is not loaded! Load the genome by calling 'library(", reference.genome, ")'")
		}else{
			genome <- get(ls(paste("package:", reference.genome, sep="")))
		}
}

checkFilesExist <- function(paths) {
  for (f in paths)
    if (! file.exists(f)) stop("Could not locate input file ", f)
}

setMethod("getCTSS",
signature(object = "CAGEset"),
function (object, sequencingQualityThreshold = 10, mappingQualityThreshold = 20, removeFirstG = TRUE, correctSystematicG = TRUE){
	
	objName <- deparse(substitute(object))
			
	sample.labels <- sampleLabels(object)
	names(sample.labels) <- rainbow(n = length(sample.labels))

	if(inputFilesType(object) == "bam" | inputFilesType(object) == "bamPairedEnd") {
		
	  checkRefGenomeIsLoaded(genomeName(object))
		
		bam.files <- inputFiles(object)
		checkFilesExist(bam.files)
					
		library.sizes <- vector()
		
		first <- TRUE
		
    param <- ScanBamParam( what = c("rname", "strand", "pos", "seq", "qual", "mapq")
                         , flag = scanBamFlag(isUnmappedQuery = FALSE))
    if (inputFilesType(object) == "bamPairedEnd")
      bamFlag(param) <- scanBamFlag(isUnmappedQuery = FALSE, isProperPair = TRUE, isFirstMateRead = TRUE)

		for(i in 1:length(bam.files)) {
	
			message("\nReading in file: ", bam.files[i], "...")
		
			bam <- scanBam(bam.files[i], param = param)
		
			message("\t-> Filtering out low quality reads...")
		
			qa.avg <- as.integer(mean(as(bam[[1]]$qual, "IntegerList")))
		
			reads.GRanges <- GRanges(seqnames = as.vector(bam[[1]]$rname), IRanges(start = bam[[1]]$pos, width = width(bam[[1]]$seq)), strand = bam[[1]]$strand, qual = qa.avg, mapq = bam[[1]]$mapq, seq = bam[[1]]$seq, read.length = width(bam[[1]]$seq))	
			reads.GRanges <- reads.GRanges[seqnames(reads.GRanges) %in% seqnames(genome)]
			reads.GRanges <- reads.GRanges[!(end(reads.GRanges) > seqlengths(genome)[as.character(seqnames(reads.GRanges))])]
			elementMetadata(reads.GRanges)$mapq[is.na(elementMetadata(reads.GRanges)$mapq)] <- Inf
			reads.GRanges.plus <- reads.GRanges[(as.character(strand(reads.GRanges)) == "+" & elementMetadata(reads.GRanges)$qual >= sequencingQualityThreshold) & elementMetadata(reads.GRanges)$mapq >= mappingQualityThreshold]
			reads.GRanges.minus <- reads.GRanges[(as.character(strand(reads.GRanges)) == "-" & elementMetadata(reads.GRanges)$qual >= sequencingQualityThreshold) & elementMetadata(reads.GRanges)$mapq >= mappingQualityThreshold]
		
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
	
	
    }else if(inputFilesType(object) == "bed") {
    
        checkRefGenomeIsLoaded(genomeName(object))

        bed.files <- inputFiles(object)
        
        checkFilesExist(bed.files)
        
        library.sizes <- vector()
        first <- TRUE
        
        for(i in 1:length(bed.files)) {
            
            message("\nReading in file: ", bed.files[i], "...")
            
            reads.GRanges <- import.bed(con = bed.files[i])
            values(reads.GRanges) <- NULL
            
            reads.GRanges <- reads.GRanges[seqnames(reads.GRanges) %in% seqnames(genome)]
            reads.GRanges <- reads.GRanges[!(end(reads.GRanges) > seqlengths(genome)[as.character(seqnames(reads.GRanges))])]
            
            reads.GRanges.plus <- reads.GRanges[strand(reads.GRanges) == "+"]
            reads.GRanges.minus <- reads.GRanges[strand(reads.GRanges) == "-"]
            
            message("\t-> Making CTSSs and counting number of tags...")

            CTSS.plus <- data.frame(chr = as.character(seqnames(reads.GRanges.plus)), pos = as.integer(start(reads.GRanges.plus)), strand = rep("+", times = length(reads.GRanges.plus)), stringsAsFactors = F)
            CTSS.minus <- data.frame(chr = as.character(seqnames(reads.GRanges.minus)), pos = as.integer(end(reads.GRanges.minus)), strand = rep("-", times = length(reads.GRanges.minus)), stringsAsFactors = F)
            CTSS <- rbind(CTSS.plus, CTSS.minus)
            CTSS$tag_count <- 1
            CTSS <- data.table(CTSS)
            CTSS <- CTSS[, as.integer(sum(tag_count)), by = list(chr, pos, strand)]
            
            setnames(CTSS, c("chr", "pos", "strand", sample.labels[i]))
            setkey(CTSS, chr, pos, strand)
        
            library.sizes <- c(library.sizes, as.integer(sum(data.frame(CTSS)[,4])))
        
            if(first == TRUE) {
                CTSS.all.samples <- CTSS
            }else{
                CTSS.all.samples <- merge(CTSS.all.samples, CTSS, all.x = TRUE, all.y = TRUE)
            }
        
            first <- FALSE

        }
        

    }else if(inputFilesType(object) == "ctss") {
	
		first <- TRUE

		ctss.files <- inputFiles(object)
		
		checkFilesExist(ctss.files)
				
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
		
	}else if(inputFilesType(object) == "CTSStable"){
	
		ctss.table.file <- inputFiles(object)
		
		if(length(ctss.table.file) > 1){
			stop("Only one file should be provided when inputFilesType = \"CTSStable\"!")
		}
		
		checkFilesExist(ctss.table.file)
		
		CTSS.all.samples <- read.table(file = ctss.table.file, header = F, stringsAsFactors = FALSE)
		if(ncol(CTSS.all.samples) != (length(sample.labels) + 3)){
			stop("Number of provided sample labels must match the number of samples in the CTSS table!")
		}
		library.sizes <- as.integer(apply(CTSS.all.samples[,c(4:ncol(CTSS.all.samples)),drop=F], 2, sum))

	}else{
		
		stop("'inputFilesType' must be one of the supported file types (\"bam\", \"bamPairedEnd\", \"ctss\", \"CTSStable\")")
		
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

#' coerceInBSgenome
#' 
#' A private (non-exported) function to discard any range that is
#' not compatible with the CAGEr object's BSgenome.
#' 
#' @param gr The genomic ranges to coerce.
#' @param genome The name of a BSgenome package, which must me installed,
#'   or \code{NULL} to skip coercion.
#' 
#' @return A GRanges object in which every range is guaranteed to be compatible
#' with the given BSgenome object.  The sequnames of the GRanges are also set
#' accordingly to the BSgenome.
#' 
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomeInfoDb seqnames
#' @importFrom S4Vectors %in%

coerceInBSgenome <- function(gr, genome) {
  if (is.null(genome)) return(gr)
  requireNamespace(genome)
  genome <- getExportedValue(genome, genome)
  gr <- gr[seqnames(gr) %in% seqnames(genome)]
  gr <- gr[! end(gr) > seqlengths(genome)[as.character(seqnames(gr))]]
  seqlevels(gr) <- seqlevels(genome)
  seqinfo(gr) <- seqinfo(genome)
  gr
}

#' loadFileIntoGRanges
#' 
#' A private (non-exported) function to load from each file format supported by CAGEr
#' 
#' @param filepath The path to the file to load.
#' @param filetype The type of the file
#' 
#' @return A GRanges object, where each range represents a single nucleotide,
#' and where the score represents the number of CAGE tags starting on that
#' nucleotide.
#' 
#' @seealso import.CTSS

loadFileIntoGRanges <- function( filepath
                               , filetype = c("bam", "bed", "bedctss", "CAGEscanMolecule", "ctss")) {
  if (missing(filetype)) stop("Please specify the file type.")
  filetype <- match.arg(filetype)
  switch( filetype
        , bam              = stop("BAM format not supported yet.")
        , bed              = import.bedmolecule(filepath)
        , bedctss          = import.bedCTSS(filepath)
        , CAGEscanMolecule = import.CAGEscanMolecule(filepath)
        , ctss             = import.CTSS(filepath))
}

#' import.bedmolecule
#' 
#' Imports a BED file where each line counts for one molecule in a
#' GRanges object where each line represents one nucleotide.
#' 
#' @param filepath The path to the BED file.
#' 
#' @seealso loadFileIntoGRanges
#' 
#' @importFrom rtracklayer import.bed
#' @importFrom S4Vectors Rle
#' 
#' @examples
#' # TODO: add exmaple file
#' # import.BED(system.file("extdata", "example.bed", package = "CAGEr"))

import.bedmolecule <- function(filepath) {
  gr <- rtracklayer::import.bed(filepath)
  tb <- table(promoters(gr, 0, 1))
  gr <- as(names(tb), "GRanges")
  score(gr) <- Rle(as.integer(tb))
  gr
}

#' import.bedCTSS
#' 
#' Imports a BED file where each line represents a single base, with a score
#' counting the number of CAGE transcription start sites (CTSS).
#' 
#' @return A GRanges object where each line represents one nucleotide.
#' 
#' @param filepath The path to the BED file.
#' 
#' @seealso loadFileIntoGRanges
#' 
#' @importFrom rtracklayer import.bed
#' 
#' @examples
#' # TODO: add exmaple file
#' # import.BED(system.file("extdata", "example.bed", package = "CAGEr"))

import.bedCTSS <- function(filepath) {
  gr <- rtracklayer::import.bed(filepath)
  if (length(unique(gr)) != length(gr))
    stop("Input must not contain the same range multiple times.")
  if (any(countOverlaps(gr) != 1))
    stop("Input must not contain overlapping ranges.")
  seqlevels(gr) <- sortSeqlevels(seqlevels(gr))
  mcols(gr) <- DataFrame(score = Rle(as.integer(score(gr))))
  gr
}

#' import.CTSS
#' 
#' Imports a "CTSS" file in a GRanges object
#' 
#' @param filepath The path to the "CTSS" file.
#' 
#' Note that the format of the "CTSS" files handled in this function is not
#' the same as the FANTOM5 "CTSS" files (which are plain BED).
#' 
#' @seealso loadFileIntoGRanges
#' 
#' @examples
#' import.CTSS(system.file("extdata", "Zf.high.chr17.ctss", package = "CAGEr"))

import.CTSS <- function(filepath) {
  CTSS <- read.table( file = filepath
                      , header = F
                      , sep = "\t"
                      , col.names  = c("chr",       "pos",     "strand",    "score")
                      , colClasses = c("character", "integer", "character", "integer"))
  GRanges( seqnames = CTSS$chr
           , ranges   = IRanges(CTSS$pos, width = 1)
           , strand   = CTSS$strand
           , score    = CTSS$score)
}

#' parseCAGEscanBlocksToGrangeTSS
#' 
#' Parse a string describing a block in a CAGEscan molecule, as output by
#' the "CAGEscan 3.0" pipeline.
#' 
#' @param block A character string representing a block in a CAGEscan molecule.
#' 
#' @return A GRanges object representing a TSS.
#' 
#' In CAGEscan molecules, blocks are separated by \sQuote{|}, \sQuote{,} or
#' \sQuote{;} for gap of coverage, splice junction (confident) and splice
#' junction (maybe) respectively.  Strand is "+" if first coordinate is lower
#' than the second one, and "-" otherwise.
#' 
#' @seealso import.CAGEscanMolecule
#' 
#' @examples
#' myMolecule <- paste0( "chr11:66268633-66268693,"
#'                     , "chr11:66271796-66271869;"
#'                     , "chr11:66272156-66272252|"
#'                     , "chr11:66272364-66272460")
#' myFirstBlock <- sub("[,;|].*", "", myMolecule)
#' 
#' parseCAGEscanBlocksToGrangeTSS(myFirstBlock)

parseCAGEscanBlocksToGrangeTSS <- function (blocks) {
  blocks <- strsplit(blocks, "[:-]")
  chr    <- unlist(lapply(blocks, `[[`, 1))
  fst    <- as.integer(unlist(lapply(blocks, `[[`, 2)))
  snd    <- as.integer(unlist(lapply(blocks, `[[`, 3)))
  strand <- ifelse(fst < snd, "+", "-")
  start  <- pmin(fst, snd)
  GRanges(chr, IRanges(start, w = 1), strand)
}

#' import.CAGEscanMolecule
#' 
#' Imports a CAGEscan \dQuote{molecule} file in a GRanges object
#' 
#' @param filepath The path to the \dQuote{molecule} file.
#' 
#' @seealso parseCAGEscanBlocksToGrangeTSS
#' 
#' @examples
#' # TODO import.CAGEscanMolecule(system.file("extdata", "example.molecule.txt", package = "CAGEr"))
#' 
#' @importFrom data.table fread

import.CAGEscanMolecule <- function(filepath) {
  molecules <- unname(unlist(data.table::fread(select = 9, paste( "grep -v \\#", filepath))))
  molecules <- sub("[,;|].*", "", molecules)
  parseCAGEscanBlocksToGrangeTSS(molecules)
}

setMethod( "getCTSS"
         , signature(object = "CAGEexp")
         , function ( object
                    , sequencingQualityThreshold = 10
                    , mappingQualityThreshold = 20
                    , removeFirstG = TRUE
                    , correctSystematicG = TRUE ){

  objName <- deparse(substitute(object))
  sample.labels <- sampleLabels(object)
  ctss.files <- inputFiles(object)
  
  # Step 0: Test existance of each file before spending time loading them.
  
  checkFilesExist(ctss.files)
  
  # Step 1: Load each file as GRangesList where each GRange is a CTSS data.
  
  l <- GRangesList()
  
  for (i in seq_along(ctss.files)) {
    message("\nReading in file: ", ctss.files[i], "...")
    gr <- loadFileIntoGRanges(ctss.files[i], inputFilesType(object)[i])
  gr <- coerceInBSgenome(gr, genomeName(object))
    l[[i]] <- gr
  }
  
  # Step 2: Create GRanges representing all the nucleotides with CAGE counts in the list.

  rowRanges <- sort(unique(unlist(l)))
  mcols(rowRanges) <- NULL

  # Step 3: Fold the GRangesList in a expression DataFrame of Rle-encoded counts.
  
  assay <- DataFrame(V1 = Rle(rep(0L, length(rowRanges))))
  
  expandRange <- function(global, local) {
    x <- Rle(rep(0L, length(global)))
    x[global %in% local] <- score(local)
    x
  }
  
  for (i in seq_along(l))
    assay[,i] <- expandRange(rowRanges, l[[i]])
  
  colnames(assay) <- sampleLabels(object)
  
  # Setp 4: Put the data in the appropriate slot of the MultiAssayExperiment.

  CTSStagCountSE(object) <-
    SummarizedExperiment( rowRanges = rowRanges
                        , assay = SimpleList(counts = assay))
  
  # Step 5: update the sample metadata (colData).
  
  librarySizes(object) <- unlist(lapply(CTSStagCountDF(object), sum))
  
  # Setp 6: overwrite the object in the parent environment.
  
  assign(objName, object, envir = parent.frame())
  invisible(1)
})

#' importPublicData
#' @noRd
#' @export

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


#' setColors
#' @noRd
#' @export

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


