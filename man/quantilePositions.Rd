\name{quantilePositions}
\alias{quantilePositions}
\alias{quantilePositions,CAGEset-method}
\title{Determining positions of CAGE signal quantiles within genomic region}
\description{
Calculates positions of quantiles of CAGE signal along tag clusters or consensus clusters in each CAGE dataset within CAGEset object.  The function calculates positions of both "lower" and "upper" quantiles as described in Details.
}
\usage{
quantilePositions(object, clusters, qLow = 0.1, qUp = 0.9, 
                  useMulticore = FALSE, nrCores = NULL)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
  \item{clusters}{
	Which clusters should be used.  Can be either \code{"tagClusters"} to calculate positions of quantiles in tag clusters (different set of genomic coordinates for every CAGE experiment) or \code{"consensusClusters"} to calculate positions of quantiles in consensus clusters (same set of genomic coordinates for every CAGE experiment).
}
  \item{qLow}{
	Which "lower" quantiles should be calculated.  It has to be a numeric vector of values in range [0,1].  See Details.
}
  \item{qUp}{
	Which "upper" quantiles should be calculated.  It has to be a numeric vector of values in range [0,1].  See Details.
}
  \item{useMulticore}{
	Logical, should multicore be used.  \code{useMulticore = TRUE} is supported only on Unix-like platforms.
}
  \item{nrCores}{
	Number of cores to use when \code{useMulticore = TRUE}.  Default value \code{NULL} uses all detected cores.
}
}
\details{
Position of the "lower" quantile \code{qLow} is defined as a point that divides the genomic region into two parts, so that the 5' part contains \code{< qLow * 100\%} of the CAGE signal of that region. Accordingly, position of the "upper" quantile \code{qUp} is defined as a point that divides the genomic region into two parts so that the 5' part contains \code{>= qUp * 100\%} of the CAGE signal of that region. Positions of one "lower" and one "upper" quantile (when \code{qLow <= qUp}) define a central part of the genomic region that contains \code{>= (qUp - qLow) * 100\%} of the CAGE signal of that region. Width of that central part is refered to as "interquantile width", which is a more robust measure of the promoter width than the total span of the region.
}
\value{
	When \code{clusters = "tagClusters"}, the slots \code{tagClustersQuantileLow} and \code{tagClustersQuantileUp} of the provided \code{\link{CAGEset}} object will be occupied with the positions of specified quantiles in all tag clusters for all CAGE datasets. When \code{clusters = "consensusClusters"} the slots \code{consensusClustersQuantileLow} and \code{consensusClustersQuantileUp} will be occupied by the corresponding information for consensus clusters.
}
\author{Vanja Haberle}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

quantilePositions(object = exampleCAGEset, clusters = "tagClusters",
qLow = c(0.1,0.2), qUp = c(0.8,0.9))
}



