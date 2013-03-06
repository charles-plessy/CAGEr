\name{mergeSamples}
\alias{mergeSamples}
\alias{mergeSamples,CAGEset,numeric-method}
\title{Merging CAGE datasets (samples)}
\description{
Merges individual CAGE datasets (experiments, samples) within the CAGEset object into specified groups. 
}
\usage{
mergeSamples(object, mergeIndex, mergedSampleLabels)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
  \item{mergeIndex}{
	Integer vector specifying which experiments should be merged.  Must be the same length as the number of datasets in the CAGEset object.  See Details.
}
  \item{mergedSampleLabels}{
	Labels for the merged datasets.  Must be the same lenght as the number of unique values in \code{mergeIndex}.  See Details.
}
}
\details{
	This function merges CAGE datasets within the CAGEset object at the level of CTSS tag counts, \emph{i.e.} tag counts of individual CTSS are summed over a group of datasets that are being merged.  After merging, all other slots in the CAGEset object will be reset and any previous data for individual experiments will be removed.
	\cr
	\code{mergeIndex} controls which datasets will be merged.  It is an integer vector that assigns a value to each dataset in the CAGEset object in the same order as they are returned by \code{sampleLabels(object)}.  Datasets with the same value in the \code{mergeIndex} will be merged.  For example, if there are 8 CAGE datasets in the CAGEset object and \code{mergeIndex = c(1,1,2,2,3,2,4,4)}, this will merge a) samples 1 and 2 b) samples 3, 4 and 6 c) samples 7 and 8, and it will leave sample 5 as it is, resulting in 4 final merged datasets. 
	\cr
	Labels provided in \code{mergedSampleLabels} will be assigned to merged datasets in the ascending order of \code{mergeIndex} values, \emph{i.e.} first label will be assigned to a dataset created by merging datasets labeled with lowest \code{mergeIndex} value (in this case \code{1}), \emph{etc}.
}
\value{
The slots \code{sampleLabels}, \code{librarySizes} and \code{tagCountMatrix} of the provided \code{\link{CAGEset}} object will be updated with the information on merged CAGE datasets and will replace the previous information on individual CAGE datasets.  All further slots with downstream information will be reset.
}
\author{Vanja Haberle}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

mergeSamples(exampleCAGEset, mergeIndex = c(1,1,2),
mergedSampleLabels = c("mergedSample1", "mergedSample2"))
exampleCAGEset
}
