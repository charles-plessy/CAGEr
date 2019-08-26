#' @include AllClasses.R CAGEr.R CAGEexp.R ClusteringMethods.R CTSS.R GetMethods.R SetMethods.R

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
#' @param useMulticore Logical, should multicore be used.
#'   \code{useMulticore = TRUE} has no effect on non-Unix-like platforms.
#'   
#' @param nrCores Number of cores to use when \code{useMulticore = TRUE}
#'   (set to \code{NULL} to use all detected cores).
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
#' the expression data can be retreived with \code{\link{CTSStagCount}} functions.  In addition,
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
#' @importFrom  S4Vectors DataFrame
#' @importFrom  S4Vectors Rle
#' @export

setGeneric( "getCTSS"
          , function( object
                    , sequencingQualityThreshold = 10, mappingQualityThreshold = 20
                    , removeFirstG = TRUE, correctSystematicG = TRUE
                    , useMulticore = FALSE, nrCores = NULL)
	            standardGeneric("getCTSS"))

checkFilesExist <- function(paths) {
  for (f in paths)
    if (! file.exists(f)) stop("Could not locate input file ", f)
}

#' toCTSSdt
#' 
#' Non-exported private helper function
#' 
#' @return Returns a \code{data.table} with TSS indicated in the \dQuote{chr},
#' \dQuote{pos} and \dQuote{strand} columns, and a score column named with
#' the sample label.
#' 
#' @noRd
#' @importFrom data.table data.table
#' @importFrom data.table setkeyv
#' @importFrom data.table setnames

toCTSSdt <- function(reads.GRanges, sample.label) {
  reads.GRanges.plus  <- reads.GRanges[strand(reads.GRanges) == "+"]
  reads.GRanges.minus <- reads.GRanges[strand(reads.GRanges) == "-"]
  gr2ctssdf <- function(gr, strand)
    data.frame( chr    = as.character(seqnames(gr))
              , pos    = if (strand == "+") start(gr) else end(gr)
              , strand = strand
              , stringsAsFactors = F)
  CTSS.plus <- gr2ctssdf(reads.GRanges.plus, "+")
	CTSS.minus <-gr2ctssdf(reads.GRanges.minus, "-")
  CTSS <- rbind(CTSS.plus, CTSS.minus)
	CTSS$tag_count <- 1L
	CTSS <- .byCtss(data.table(CTSS), "tag_count", sum)
	setnames(CTSS, c("chr", "pos", "strand", sample.label))
	setkeyv(CTSS, c("chr", "pos", "strand"))
	CTSS
}

addCTSScolumn <- function(CTSS.all.samples, CTSS) {
  if(is.null(CTSS.all.samples))
    return(CTSS)
  merge(CTSS.all.samples, CTSS, all.x = TRUE, all.y = TRUE)
}

#' @rdname getCTSS

setMethod("getCTSS", "CAGEset", function ( object
                                         , sequencingQualityThreshold, mappingQualityThreshold
                                         , removeFirstG, correctSystematicG
                                         , useMulticore, nrCores) {
  # Initialise values
  objName <- deparse(substitute(object))
  sample.labels <- sampleLabels(object)
  names(sample.labels) <- rainbow(n = length(sample.labels))
  CTSS.all.samples <- NULL
  checkFilesExist(inputFiles(object))
  
  # Switch on file types

  if(inputFilesType(object) == "bam" | inputFilesType(object) == "bamPairedEnd") {
    genome <- getRefGenome(genomeName(object))
    for(i in 1:length(inputFiles(object))) {
      message("\nReading in file: ", inputFiles(object)[i], "...")
      reads.GRanges <- import.bam( filepath                   = inputFiles(object)[i]
                                 , filetype                   = inputFilesType(object)
                                 , sequencingQualityThreshold = sequencingQualityThreshold
                                 , mappingQualityThreshold    = mappingQualityThreshold)
      reads.GRanges <- coerceInBSgenome(reads.GRanges, genomeName(object))
      if(removeFirstG == TRUE){
        CTSS <- .remove.added.G(reads.GRanges, genome, correctSystematicG = correctSystematicG, sample.label = sample.labels[i])
      }else{
        CTSS <- toCTSSdt(reads.GRanges, sample.labels[i])
      }
      message("\t-> Making CTSSs and counting number of tags...")
      CTSS.all.samples <- addCTSScolumn(CTSS.all.samples, CTSS)
    }

  }else if(inputFilesType(object) == "bed") {
    
    genome <- getRefGenome(genomeName(object))
    for(i in 1:length(inputFiles(object))) {
      message("\nReading in file: ", inputFiles(object)[i], "...")
      reads.GRanges <- import.bed(con = inputFiles(object)[i])
      values(reads.GRanges) <- NULL
      reads.GRanges <- coerceInBSgenome(reads.GRanges, genomeName(object))
      message("\t-> Making CTSSs and counting number of tags...")
      CTSS <- toCTSSdt(reads.GRanges, sample.labels[i])
      CTSS.all.samples <- addCTSScolumn(CTSS.all.samples, CTSS)
    }

  }else if(inputFilesType(object) == "ctss") {
    
    for(i in 1:length(inputFiles(object))) {
      message("\nReading in file: ", inputFiles(object)[i], "...")
      CTSS <- read.table(file = inputFiles(object)[i], header = F, sep = "\t", colClasses = c("character", "integer", "character", "integer"), col.names = c("chr", "pos", "strand", sample.labels[i]))
      CTSS <- data.table(CTSS)
      setkeyv(CTSS, cols = c("chr", "pos", "strand"))
      CTSS.all.samples <- addCTSScolumn(CTSS.all.samples, CTSS)
    }
  
  }else if(inputFilesType(object) == "CTSStable"){
    
    if(length(inputFiles(object)) > 1)
      stop("Only one file should be provided when inputFilesType = \"CTSStable\"!")
    CTSS.all.samples <- read.table(file = inputFiles(object), header = F, stringsAsFactors = FALSE)
    if(ncol(CTSS.all.samples) != (length(sample.labels) + 3))
      stop("Number of provided sample labels must match the number of samples in the CTSS table!")

  }else{
	  
    stop("'inputFilesType' must be one of the supported file types (\"bam\", \"bamPairedEnd\", \"ctss\", \"CTSStable\")")
  }
  
  # Brush up the CTSS table
  
  CTSS.all.samples <- data.frame(CTSS.all.samples)
  for(i in 4:ncol(CTSS.all.samples)){
    CTSS.all.samples[is.na(CTSS.all.samples[,i]),i] <- as.integer(0)
  }
  CTSS.all.samples <- CTSS.all.samples[order(CTSS.all.samples$chr, CTSS.all.samples$pos),]
  rownames(CTSS.all.samples) <- c(1:nrow(CTSS.all.samples))
	
	# Update the object
	
  object@sampleLabels <- sample.labels
  object@CTSScoordinates <- CTSS.all.samples[,c("chr", "pos", "strand")]
  object@tagCountMatrix <- as.data.frame(CTSS.all.samples[,c(4:ncol(CTSS.all.samples)),drop=F])
  cat("\n")
  object@librarySizes <- as.integer(colSums(object@tagCountMatrix))
  names(object@librarySizes) <- object@sampleLabels
  assign(objName, object, envir = parent.frame())
  invisible(1)
})

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
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomeInfoDb seqlevels seqlevels<-
#' @importFrom GenomeInfoDb seqnames
#' @importFrom S4Vectors %in%

coerceInBSgenome <- function(gr, genome) {
  if (is.null(genome)) return(gr)
  genome <- getRefGenome(genome)
  gr <- gr[seqnames(gr) %in% seqnames(genome)]
  gr <- gr[! end(gr) > seqlengths(genome)[as.character(seqnames(gr))]]
  seqlevels(gr) <- seqlevels(genome)
  seqinfo(gr) <- seqinfo(genome)
  gr
}

#' loadFileIntoGPos
#' 
#' A private (non-exported) function to load from each file format supported by CAGEr
#' 
#' @param filepath The path to the file to load.
#' @param filetype The type of the file
#' @param sequencingQualityThreshold See getCTSS().
#' @param mappingQualityThreshold See getCTSS().
#' @param removeFirstG See getCTSS().
#' @param correctSystematicG See getCTSS().
#' @param genome See coerceInBSgenome().
#' 
#' @return A [GPos()] object where the score represents the number of CAGE tags
#' starting on that nucleotide.
#' 
#' @seealso import.CTSS
#' @family loadFileIntoGPos

loadFileIntoGPos    <- function( filepath
                               , filetype = c( "bam", "bamPairedEnd"
                                             , "bed", "bedctss", "bedScore"
                                             , "CAGEscanMolecule", "ctss")
                               , sequencingQualityThreshold
                               , mappingQualityThreshold
                               , removeFirstG
                               , correctSystematicG
                               , genome) {
  if (missing(filetype)) stop("Please specify the file type.")
  filetype <- match.arg(filetype)
  switch( filetype
        , bam              = import.bam.ctss( filepath = filepath
                                            , filetype = "bam"
                                            , sequencingQualityThreshold = sequencingQualityThreshold
                                            , mappingQualityThreshold = mappingQualityThreshold
                                            , removeFirstG = removeFirstG
                                            , correctSystematicG = correctSystematicG
                                            , genome = genome)
        , bamPairedEnd     = import.bam.ctss( filepath = filepath
                                            , filetype = "bamPairedEnd"
                                            , sequencingQualityThreshold = sequencingQualityThreshold
                                            , mappingQualityThreshold = mappingQualityThreshold
                                            , removeFirstG = removeFirstG
                                            , correctSystematicG = correctSystematicG
                                            , genome = genome)
        , bed              = import.bedmolecule(filepath)
        , bedScore         = import.bedScore(filepath)
        , bedctss          = import.bedCTSS(filepath)
        , CAGEscanMolecule = import.CAGEscanMolecule(filepath)
        , ctss             = import.CTSS(filepath))
}

#' moleculesGR2CTSS
#' 
#' Calculates CTSS positions from a GenomicRanges object where each element
#' represents a single molecule.
#' 
#' @param gr A [GRanges] object.
#' 
#' @return Returns a [GRanges] object.
#' 
#' @family loadFileIntoGPos
#' 
#' @importFrom S4Vectors Rle
#' 
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10), c("+", "-", "+"))
#' CAGEr:::moleculesGR2CTSS(gr)

moleculesGR2CTSS <- function(gr) {
  tb <- table(promoters(gr, 0, 1))
  gp <- as(names(tb), "UnstitchedGPos")
  score(gp) <- Rle(as.integer(tb))
  gp
}

#' import.bam
#' 
#' Imports CTSS data from a BAM file.
#' 
#' @param filepath The path to the BAM file.
#' @param filetype bam or bamPairedEnd.
#' @param sequencingQualityThreshold See getCTSS().
#' @param mappingQualityThreshold See getCTSS().
#' 
#' @family loadFileIntoGPos
#' 
#' @importFrom Rsamtools bamFlag<-
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools scanBamFlag
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomicAlignments qwidth
#' 
#' @examples
#' # TODO: add exmaple file
#' # import.bam(system.file("extdata", "example.bam", package = "CAGEr"))

import.bam <- function( filepath
                      , filetype
                      , sequencingQualityThreshold = 10
                      , mappingQualityThreshold = 20) {
  param <- ScanBamParam( what       = c("seq", "qual")
                       , flag       = scanBamFlag(isUnmappedQuery = FALSE)
                       , mapqFilter = mappingQualityThreshold)
  if (filetype == "bamPairedEnd")
    bamFlag(param) <- scanBamFlag( isUnmappedQuery = FALSE
                                 , isProperPair    = TRUE
                                 , isFirstMateRead = TRUE)
  bam <- readGAlignments(filepath, param = param)
  message("\t-> Filtering out low quality reads...")
  # We need to loop on qualities because there is a hard limit on the length.
  # See <https://support.bioconductor.org/p/97993/>.
  qual <- mcols(bam)$qual
  start <- 1
  chunksize <- 1e6
  qual.avg <- vector(mode = "integer")
  repeat {
    if (start + chunksize <= length(qual)) {
      end <- start + chunksize
    } else {
      end <- length(qual)
    }
    qual.avg <- c(qual.avg, as.integer(mean(as(qual[start:end], "IntegerList"))))
    if (end == length(qual)) {
      break
    } else {
      start <- end + 1
    }
  }
  bam <- bam[qual.avg >= sequencingQualityThreshold]
  gr <- GRanges(bam)
  gr$read.length <- qwidth(bam)
  gr
}

#' bam2CTSS
#' 
#' Converts from BAM to CTSS
#' 
#' Converts genomic ranges representing SAM/BAM alignments into a CTSS object.
#' 
#' @param gr A \code{\link{GRanges}} object returned by \code{\link{import.bam}}.
#' @param removeFirstG See getCTSS().
#' @param correctSystematicG See getCTSS().
#' @param genome See coerceInBSgenome().
#' 
#' @return Returns a \code{\link{CTSS}} object.
#' 
#' @importFrom GenomeInfoDb bsgenomeName
#' @family loadFileIntoGPos

bam2CTSS <- function(gr, removeFirstG, correctSystematicG, genome) {
  gr <- coerceInBSgenome(gr, genome)
  if(removeFirstG == TRUE) {
    message("\t-> Removing the first base of the reads if 'G' and not aligned to the genome...")
    gr <- .remove.added.G.CTSS(gr, genome, correctSystematicG = correctSystematicG)
  }
  tb <- table(promoters(gr, 0, 1))
  gp <- as(names(tb), "UnstitchedGPos")
  score(gp) <- Rle(unclass(tb))
  gp
} 

#' import.bam.ctss
#' 
#' Imports CTSS data from a BAM file.
#' 
#' @param filepath The path to the BAM file.
#' @param filetype bam or bamPairedEnd.
#' @param sequencingQualityThreshold See getCTSS().
#' @param mappingQualityThreshold See getCTSS().
#' @param removeFirstG See getCTSS().
#' @param correctSystematicG See getCTSS().
#' @param genome See coerceInBSgenome().
#' 
#' @return Returns a [CTSS] object.
#' 
#' @family loadFileIntoGPos

import.bam.ctss <- function( filepath, filetype, sequencingQualityThreshold
                           , mappingQualityThreshold, removeFirstG
                           , correctSystematicG, genome) {
  bam <- import.bam( filepath = filepath
                   , filetype = filetype
                   , sequencingQualityThreshold = sequencingQualityThreshold
                   , mappingQualityThreshold = mappingQualityThreshold)
  bam2CTSS( gr= bam
          , removeFirstG = removeFirstG
          , correctSystematicG = correctSystematicG
          , genome = genome)
}

#' import.bedmolecule
#' 
#' Imports a BED file where each line counts for one molecule in a
#' GRanges object where each line represents one nucleotide.
#' 
#' @param filepath The path to the BED file.
#' 
#' @return Returns a \code{\link{CTSS}} object.
#' 
#' @family loadFileIntoGPos
#' 
#' @importFrom rtracklayer import.bed
#' @importFrom S4Vectors Rle
#' 
#' @examples
#' # TODO: add exmaple file
#' # import.BED(system.file("extdata", "example.bed", package = "CAGEr"))

import.bedmolecule <- function(filepath) {
  moleculesGR2CTSS(rtracklayer::import.bed(filepath))
}

#' import.bedScore
#' 
#' Imports a BED file where the score indicates a number of counts
#' for a given alignment.
#' 
#' @param filepath The path to the BED file.
#' 
#' @return A GRanges object where each line represents one nucleotide.
#' 
#' @family loadFileIntoGPos
#' 
#' @importFrom rtracklayer import.bed
#' @importFrom S4Vectors Rle
#' 
#' @examples
#' # TODO: add exmaple file
#' # import.bedScore(system.file("extdata", "example.bed", package = "CAGEr"))

import.bedScore <- function(filepath) {
  gr <- rtracklayer::import.bed(filepath)
  df <- data.frame( p = as.character(promoters(gr, 0, 1))
                  , s = score(gr))
  df <- rowsum(df$s, df$p)
  gp <- as(rownames(df), "UnstitchedGPos")
  score(gp) <- Rle(as.integer(df[,1]))
  gp
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
#' @family loadFileIntoGPos
#' 
#' @importFrom rtracklayer import.bed
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom GenomicRanges countOverlaps
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
  as(gr, "UnstitchedGPos")
}

#' import.CTSS
#' 
#' Imports a "CTSS" file in a [GPos] object
#' 
#' @param filepath The path to the "CTSS" file.
#' 
#' Note that the format of the "CTSS" files handled in this function is not
#' the same as the FANTOM5 "CTSS" files (which are plain BED).
#' 
#' @family loadFileIntoGPos
#' 
#' @importFrom utils read.table
#' 
#' @examples
#' CAGEr:::import.CTSS(system.file("extdata", "Zf.high.chr17.ctss", package = "CAGEr"))

import.CTSS <- function(filepath) {
  CTSS <- read.table( file = filepath
                      , header = F
                      , sep = "\t"
                      , col.names  = c("chr",       "pos",     "strand",    "score")
                      , colClasses = c("character", "integer", "character", "integer"))
  gp <- GPos( stitch = FALSE
            , GRanges( seqnames = CTSS$chr
                     , ranges   = IRanges(CTSS$pos, width = 1)
                     , strand   = CTSS$strand))
  score(gp) <- CTSS$score
  gp
}

#' parseCAGEscanBlocksToGrangeTSS
#' 
#' Parse a string describing a block in a CAGEscan molecule, as output by
#' the "CAGEscan 3.0" pipeline.
#' 
#' @param blocks A character string representing a block in a CAGEscan molecule.
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
#' CAGEr:::parseCAGEscanBlocksToGrangeTSS(myFirstBlock)

parseCAGEscanBlocksToGrangeTSS <- function (blocks) {
  blocks <- strsplit(blocks, "[:-]")
  chr    <- unlist(lapply(blocks, `[[`, 1))
  fst    <- as.integer(unlist(lapply(blocks, `[[`, 2)))
  snd    <- as.integer(unlist(lapply(blocks, `[[`, 3)))
  strand <- ifelse(fst < snd, "+", "-")
  start  <- pmin(fst, snd)
  GPos(stitch=FALSE, GRanges(chr, IRanges(start, width = 1), strand))
}

#' import.CAGEscanMolecule
#' 
#' Imports a CAGEscan \dQuote{molecule} file in a [`GRanges`] object
#' 
#' @param filepath The path to the \dQuote{molecule} file.
#' 
#' @seealso parseCAGEscanBlocksToGrangeTSS
#' 
#' @importFrom GenomicRanges GRangesList
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

#' @rdname getCTSS

setMethod( "getCTSS", "CAGEexp"
         , function ( object
                    , sequencingQualityThreshold, mappingQualityThreshold
                    , removeFirstG, correctSystematicG
                    , useMulticore, nrCores) {

  objName <- deparse(substitute(object))

  # Step 0: Test existance of each file before spending time loading them.
  
  checkFilesExist(inputFiles(object))
  
  # Step 1: Load each file as GRangesList where each GRange is a CTSS data.
  
  l <- as.list(seq_along(inputFiles(object)))

  l <- bplapply(l, function(i) {
      message("\nReading in file: ", inputFiles(object)[i], "...")
      gr <-    loadFileIntoGPos( filepath                   = inputFiles(object)[i]
                               , filetype                   = inputFilesType(object)[i]
                               , sequencingQualityThreshold = sequencingQualityThreshold
                               , mappingQualityThreshold    = mappingQualityThreshold
                               , removeFirstG               = removeFirstG
                               , correctSystematicG         = correctSystematicG
                               , genome                     = genomeName(object))
      gr <- coerceInBSgenome(gr, genomeName(object))
      l[[i]] <- sort(gr)
    }, BPPARAM = CAGEr_Multicore(useMulticore, nrCores))
  
  l <-  GRangesList(l)
  
  # Step 2: Create GPos representing all the nucleotides with CAGE counts in the list.

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
  
  rowRanges <- .CTSS(rowRanges, bsgenomeName = genomeName(object))
  
  colnames(assay) <- sampleLabels(object)
  
  # Setp 4: Put the data in the appropriate slot of the MultiAssayExperiment.

  CTSStagCountSE(object) <-
    SummarizedExperiment( rowRanges = rowRanges
                        , assays    = SimpleList(counts = assay))
  
  # Step 5: update the sample metadata (colData).
  
  librarySizes(object) <- unlist(lapply(CTSStagCountDF(object), sum))
  
  # Setp 6: overwrite the object in the parent environment.
  
  assign(objName, object, envir = parent.frame())
  invisible(1)
})

#' importPublicData
#' @noRd
#' @importFrom utils data
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
			ENCODEtissueCAGEfly <- NULL
			data("ENCODEtissueCAGEfly", envir = environment())
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
			ENCODEhumanCellLinesSamples <- NULL
			data("ENCODEhumanCellLinesSamples", envir = environment())
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
		FANTOMhumanSamples <- FANTOMmouseSamples <- NULL
		data("FANTOMhumanSamples", envir = environment())
		info.df1 <- FANTOMhumanSamples
		data("FANTOMmouseSamples", envir = environment())
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
		  FANTOM5humanSamples <- NULL
			data("FANTOM5humanSamples", envir = environment())
			samples.info <- FANTOM5humanSamples
			genome.name <- "BSgenome.Hsapiens.UCSC.hg19"
		}else if(dataset == "mouse"){
		  FANTOM5mouseSamples <- NULL
			data("FANTOM5mouseSamples", envir = environment())
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
		
	  ZebrafishSamples <- NULL
		data("ZebrafishSamples", envir = environment())
		if(dataset == "ZebrafishCAGE"){
			if(group == "development"){
				if(!(all(sample %in% ZebrafishSamples$sample))){
					stop("Some sample names cannot be found for the specified dataset! Call data(ZebrafishSamples) and check the 'sample' column for valid sample names!")
				}else{
					genome.name <- "BSgenome.Drerio.UCSC.danRer7"
					ZebrafishCAGE <- NULL
					data("ZebrafishCAGE", envir = environment())
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