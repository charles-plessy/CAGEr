#' @include SetMethods.R
#' 
#' @name mergeSamples
#' 
#' @title Merge CAGE samples
#' 
#' @description Merges individual CAGE samples (datasets, experiments) within
#' the CAGEr object into specified groups.
#'  
#' @param object A \code{\link{CAGEr}} object.
#' 
#' @param mergeIndex Integer vector specifying which experiments should be merged.
#'        (one value per sample, see Details).
#'        
#' @param mergedSampleLabels Labels for the merged datasets (same length as the
#'        number of unique values in \code{mergeIndex})
#' 
#' @details The samples within the CAGEr object are merged by adding the raw tag
#' counts of individual CTSS that belong tho the same group.  After merging, all
#' other slots in the CAGEr object will be reset and any previous data for
#' individual experiments will be removed.
#' 
#' \code{mergeIndex} controls which samples will be merged.  It is an integer
#' vector that assigns a group identifier to each sample, in the same order as
#' they are returned by \code{sampleLabels(object)}.  For example, if there are
#' 8 CAGE samples in the CAGEr object and \code{mergeIndex = c(1,1,2,2,3,2,4,4)},
#' this will merge a) samples 1 and 2, b) samples 3, 4 and 6,  c) samples 7 and
#' 8, and d) it will leave sample 5 as it is, resulting in 4 final merged datasets.
#' 
#' Labels provided in \code{mergedSampleLabels} will be assigned to merged datasets in the ascending
#' order of \code{mergeIndex} values, \emph{i.e.} first label will be assigned to a dataset created
#' by merging datasets labeled with lowest \code{mergeIndex} value (in this case \code{1}),
#' \emph{etc}.
#' 
#' @return The slots \code{sampleLabels}, \code{librarySizes} and \code{tagCountMatrix} of the
#' provided \code{\link{CAGEr}} object will be updated with the information on merged CAGE datasets
#' and will replace the previous information on individual CAGE datasets.  All further slots with
#' downstream information will be reset.
#' 
#' @author Vanja Haberle
#' @author Charles Plessy
#' 
#' @examples 
#' 
#' mergeSamples( exampleCAGEexp
#'             , mergeIndex = c(3,2,4,4,1)
#'             , mergedSampleLabels = c("zf_unfertilized", "zf_high", "zf_30p_dome", "zf_prim6"))
#' exampleCAGEexp
#' 
#' @export

setGeneric("mergeSamples", function(object, mergeIndex, mergedSampleLabels)
	standardGeneric("mergeSamples"))

checkMergeOK <- function(object, mergeIndex, mergedSampleLabels) {
  if( length(mergeIndex) != length(sampleLabels(object)))
    stop( "length of ", sQuote("mergeIndex"), " must match number of samples! See the "
        ,  sQuote("sampleLabels()"), "function to list your CAGE samples.")

  if( length(unique(mergeIndex)) != length(mergedSampleLabels))
    stop( "numer of provided ", sQuote("mergedSampleLabels"), " must match number of unique values "
        , "provided in ", sQuote("mergeIndex"), "!")

  if(! identical(mergedSampleLabels, make.names(mergedSampleLabels)))
    stop("'mergedSampleLabels' must contain non-empty strings beginning with a letter!")
	
  if(length(unique(mergedSampleLabels)) != length(mergedSampleLabels))
    stop("Duplicated sample labels are not allowed!")
}

#' @rdname mergeSamples

setMethod( "mergeSamples", "CAGEexp", function (object, mergeIndex, mergedSampleLabels) {
  checkMergeOK(object, mergeIndex, mergedSampleLabels)

  tag.count.DF.new <- tapply(
    as.list(CTSStagCountDF(object)),
    mergeIndex,
    \(l) Reduce( f    = `+`
               , x    = DataFrame(l)
               , init = l |> DataFrame() |> nrow() |> rep(x=0L) |> Rle())
  )
  tag.count.DF.new <- DataFrame(do.call(list, tag.count.DF.new))
  colnames(tag.count.DF.new) <- mergedSampleLabels

  # Use c() to force tapply to return a plain vector (not an array)
  lib.sizes.new <- c(tapply(librarySizes(object), mergeIndex, sum))
  inputFilesType.new <- c(tapply(inputFilesType(object), mergeIndex, function(X) paste(unique(X))))
  
  new.CAGE.exp <- DataFrame( inputFiles     = paste(mergedSampleLabels, "_merged", sep = "")
                           , inputFilesType = inputFilesType.new
                           , sampleLabels   = mergedSampleLabels
                           , librarySizes   = lib.sizes.new)
  rownames(new.CAGE.exp) <- mergedSampleLabels
  
  new.CAGE.exp <- CAGEexp( metadata = list(genomeName = genomeName(object))
                     , colData  = new.CAGE.exp)
  
  CTSStagCountSE(new.CAGE.exp) <-
    SummarizedExperiment( rowRanges = rowRanges(CTSStagCountSE(object))
                        , assays    = SimpleList(counts = tag.count.DF.new))
  
  setColors(new.CAGE.exp, rainbow(n = length(mergedSampleLabels)))

  new.CAGE.exp
})
