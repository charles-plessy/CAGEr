\name{cumulativeCTSSdistribution}
\alias{cumulativeCTSSdistribution}
\alias{cumulativeCTSSdistribution,CAGEset-method}
\title{Calculating cumulative sum of CAGE signal along genomic region}
\description{
	Calculates cumulative sum of CAGE signal along each tag cluster or consensus cluster in every CAGE dataset within CAGEset object.
}
\usage{
cumulativeCTSSdistribution(object, clusters, useMulticore = FALSE, 
                           nrCores = NULL)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
  \item{clusters}{
	Which clusters should be used.  Can be either \code{clusters = "tagClusters"} to calculate cumulative sum along tag clusters (different set of genomic coordinates for every CAGE experiment) or \code{clusters = "consensusClusters"} to calculate cumulative sum along consensus clusters (same set of genomic coordinates for every CAGE experiment).
}
  \item{useMulticore}{
	Logical, should multicore be used.  \code{useMulticore = TRUE} is supported only on Unix-like platforms.
}
  \item{nrCores}{
	Number of cores to use when \code{useMulticore = TRUE}.  Default value \code{NULL} uses all detected cores.
}
}
\value{
	The slot \code{CTSScumulativesTagClusters} (when \code{clusters = "tagClusters"}) or \code{CTSScumulativesConsensusClusters} (when \code{clusters = "consensusClusters"}) of the provided \code{\link{CAGEset}} object will be occupied by the list containing cumulative sum of the CAGE signal along genomic regions per CAGE experiment.
}
\author{Vanja Haberle}

\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

cumulativeCTSSdistribution(object = exampleCAGEset, clusters = "tagClusters")
}
