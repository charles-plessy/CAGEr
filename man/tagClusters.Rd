\name{tagClusters}
\alias{tagClusters}
\alias{tagClusters,CAGEset-method}
\title{Extracting tag clusters (TCs) for individual CAGE experiment from CAGEset object}
\description{
	Extracts tag clusters (TCs) produced by \code{\link{clusterCTSS}} function for a specified CAGE experiment from a CAGEset object.
}
\usage{
tagClusters(object, sample)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
  \item{sample}{
	Label of the CAGE dataset (experiment, sample) for which to extract tag clusters.
}
}
\value{
	Returns a \code{data.frame} with genomic coordinates, position of dominant TSS, total CAGE signal and additional information for all TCs from specified CAGE dataset.
}
\author{Vanja Haberle}
\seealso{
\code{\link{clusterCTSS}}
}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

TC <- tagClusters(object = exampleCAGEset, sample = "sample2")
head(TC)
}
