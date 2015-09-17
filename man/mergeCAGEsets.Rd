\name{mergeCAGEsets}
\alias{mergeCAGEsets}
\alias{mergeCAGEsets,CAGEset,CAGEset-method}
\title{Merging two CAGEset objects into one}
\description{
	Merges two \code{\link{CAGEset}} objects into one by combining the CTSS genomic coordinates and raw tag counts. Merged \code{\link{CAGEset}} object will contain a union of TSS positions present in the two input \code{\link{CAGEset}} objects and raw tag counts for those TSSs in all samples from both input \code{\link{CAGEset}} objects.
	}
\usage{
mergeCAGEsets(cs1, cs2)
}
\arguments{
  \item{cs1}{
	A \code{\link{CAGEset}} object
	}
  \item{cs2}{
	A \code{\link{CAGEset}} object
	}
}
\details{
Note that merging discards all other information present in the two \code{\link{CAGEset}} objects, \emph{i.e.} the merged object will not contain any normalised tag counts, CTSS clusters, quantile positions, \emph{etc.}, so these have to be calculated again by calling the appropriate functions on the merged \code{\link{CAGEset}} object. Also, it is only possible to merge two \code{\link{CAGEset}} objects that contain TSS information for the same reference genome and do not share any sample names.
}
\value{
	A \code{\link{CAGEset}} object is returned, which contains a union of TSS positions present in the two input \code{\link{CAGEset}} objects and raw tag counts for those TSSs in all samples from both input \code{\link{CAGEset}} objects.
}
\author{Vanja Haberle}

\seealso{
\code{\link{CAGEset}}
}
\examples{
library(BSgenome.Drerio.UCSC.danRer7)

pathsToInputFiles <- system.file("extdata", c("Zf.unfertilized.egg.chr17.ctss",
"Zf.30p.dome.chr17.ctss", "Zf.prim6.rep1.chr17.ctss"), package="CAGEr")

myCAGEset1 <- new("CAGEset", genomeName = "BSgenome.Drerio.UCSC.danRer7",
inputFiles = pathsToInputFiles[1:2], inputFilesType = "ctss", sampleLabels = 
c("sample1", "sample2"))
getCTSS(myCAGEset1)

myCAGEset2 <- new("CAGEset", genomeName = "BSgenome.Drerio.UCSC.danRer7",
inputFiles = pathsToInputFiles[3], inputFilesType = "ctss", sampleLabels = 
"sample3")
getCTSS(myCAGEset2)

myCAGEset <- mergeCAGEsets(myCAGEset1, myCAGEset2)
}
