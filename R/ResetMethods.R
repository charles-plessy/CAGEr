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
#' @return Returns a `CAGEexp` object, which contains a union of
#' TSS positions present in the two input objects and raw tag counts for those
#' TSSs in all samples from both input objects.
#' 
#' @author Vanja Haberle
#' @author Charles Plessy
#' 
#' @seealso [`CAGEexp`]
#' 
#' @examples 
#' library(BSgenome.Drerio.UCSC.danRer7)
#' 
#' pathsToInputFiles <- system.file("extdata", c("Zf.unfertilized.egg.chr17.ctss",
#'   "Zf.30p.dome.chr17.ctss", "Zf.prim6.rep1.chr17.ctss"), package="CAGEr")
#'   
#' ce1 <- CAGEexp(genomeName = "BSgenome.Drerio.UCSC.danRer7",
#' inputFiles = pathsToInputFiles[1:2], inputFilesType = "ctss", sampleLabels =
#' c("sample1", "sample2"))
#' ce1 <- getCTSS(ce1)
#' 
#' ce2 <- CAGEexp(genomeName = "BSgenome.Drerio.UCSC.danRer7",
#' inputFiles = pathsToInputFiles[3], inputFilesType = "ctss", sampleLabels =
#' "sample3")
#' 
#' ce2 <- getCTSS(ce2)
#' 
#' ce <- mergeCAGEsets(ce1, ce2)
#' 
#' @export

setGeneric( "mergeCAGEsets"
          , function(cs1, cs2) {
            
  if (genomeName(cs1) != genomeName(cs2))
    stop("Cannot merge two objects with data from different genomes!")
  
  if(any(sampleLabels(cs1) %in% sampleLabels(cs2)))
    stop("Cannot merge two objects that share same sample labels!")
  
  standardGeneric("mergeCAGEsets")
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
