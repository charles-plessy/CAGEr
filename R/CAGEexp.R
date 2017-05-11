#' @rdname CAGEexp
#' 
#' @examples 
#' 
#' pathsToInputFiles <- list.files( system.file("extdata", package = "CAGEr")
#'                                , "ctss$"
#'                                , full.names = TRUE)
#' 
#' myCAGEexp <-
#'   new( "CAGEexp"
#'      , colData = DataFrame( inputFiles = pathsToInputFiles
#'                           , sampleLabels = sub(".chr17.ctss"
#'                                               , ""
#'                                               , basename(pathsToInputFiles)))
#'      , metadata = list( genomeName = "BSgenome.Drerio.UCSC.danRer7"
#'                       , inputFilesType = "ctss"))
#' 
#' getCTSS(myCAGEexp)
#' 
#' colData(myCAGEexp)
#' 
#' @export

setClass("CAGEexp",
  contains = "MultiAssayExperiment",
  validity = function(object) {
    #		if(!(object@genomeName %in% suppressWarnings(suppressMessages(BSgenome::available.genomes()))))
    #		return("'genomeName' must be a name of one of the genome packages available in BSgenome! See 'available.genomes()'")
    #		if(object@genomeName %in% rownames(installed.packages()) == FALSE)
    #		return("Requested genome is not installed! Please install required BSgenome package before running CAGEr.")
    
    if (is.null(metadata(object)$genomeName))
      return("Missing BSgenome.")
    
    if (is.null(colData(object)$inputFiles))
      return("Missing input file list.")
    
    if (is.null(metadata(object)$inputFilesType))
      return("Missing input file type.")
    
    if (! metadata(object)$inputFilesType %in%
          c( "bam", "bamPairedEnd", "bed", "ctss", "CTSStable"
           , "FANTOM5", "ENCODE", "FANTOM3and4", "ZebrafishDevelopment"))
      return("'inputFilesType' must be one of supported input file types (\"bam\", \"bamPairedEnd\", \"bed\", \"ctss\", \"CTSStable\")!")
    
    if (is.null(colData(object)$sampleLabels))
      return("Missing sample labels.")
    
    if (!(all(nzchar(colData(object)$sampleLabels))) |
        !(all(substr(colData(object)$sampleLabels, start = 1, stop = 1) %in%c(letters, LETTERS))))
      return("All sample labels must be a non-empty strings beginning with a letter!")
    
    if (length(unique(colData(object)$sampleLabels)) !=
        length(colData(object)$sampleLabels)) 
      return("Duplicated sample labels are not allowed!")
  }
)

#' loadFileIntoGRanges
#' 
#' A private (non-exported) function to load from each file format supported by CAGEr
#' 
#' @param filepatch The path to the file to load.
#' @param filetype The type of the file
#' 
#' @return A GRanges object, where each range represents a single nucleotide,
#' and where the score represents the number of CAGE tags starting on that
#' nucleotide.
#' 
#' @seealso import.CTSS

loadFileIntoGRanges <- function( filepath
                               , filetype = c("bam", "bed", "ctss")) {
  if (missing(filetype)) stop("Please specify the file type.")
  filetype <- match.arg(filetype)
  switch( filetype
        , bam  = stop("BAM format not supported yet.")
        , bed  = stop("BED format not supported yet.")
        , ctss = import.CTSS(filepath))
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

#' getCTSS

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
  
  for (f in ctss.files)
    if (! file.exists(f)) stop("Could not locate input file ", f)
  
  # Step 1: Load each file as GRangesList where each GRange is a CTSS data.
  
  l <- GRangesList()
  
  for (i in seq_along(ctss.files)) {
    message("\nReading in file: ", ctss.files[i], "...")
    l[[i]] <- import.CTSS(ctss.files[i])
  }
  
  # Step 2: Create GRanges representing all the nucleotides with CAGE counts in the list.

  rowRanges <- unique(unlist(l))
  rowRanges <- rowRanges[order(as.character(rowRanges))]
  
  # Step 3: Fold the GRangesList in a expression DataFrame of Rle-encoded counts.
  
  assay <- DataFrame(V1 = Rle(rep(0, length(rowRanges))))
  
  expandRange <- function(global, local) {
    x <- Rle(rep(0, length(global)))
    x[global %in% local] <- score(local)
    x
  }
  
  for (i in seq_along(l))
    assay[,i] <- expandRange(rowRanges, l[[i]])
  
  colnames(assay) <- sampleLabels(object)
  
  # Setp 4: Put the data in the appropriate slot of the MultiAssayExperiment.

  object@ExperimentList$tagCountMatrix <-
    SummarizedExperiment( rowRanges = rowRanges
                        , assay = SimpleList(counts = assay))
  
  # Step 5: update the sample metadata (colData).
  
  colData(object)$librarySizes <- unlist(lapply(assay(object@ExperimentList$tagCountMatrix), sum))
  
  # Setp 6: overwrite the object in the parent environment.
  
  assign(objName, object, envir = parent.frame())
  invisible(1)
})