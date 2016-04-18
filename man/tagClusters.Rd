\name{tagClusters}
\alias{tagClusters}
\alias{tagClusters,CAGEset-method}
\title{Extracting tag clusters (TCs) for individual CAGE experiment from CAGEset object}
\description{
	Extracts tag clusters (TCs) produced by \code{\link{clusterCTSS}} function for a specified CAGE experiment from a CAGEset object.
}
\usage{
tagClusters(object, sample, returnInterquantileWidth = FALSE, 
			qLow = NULL, qUp = NULL)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
  \item{sample}{
	Label of the CAGE dataset (experiment, sample) for which to extract tag clusters.
}
  \item{returnInterquantileWidth}{
	Should the interquantile width for each tag cluster be returned.
}
  \item{qLow}{
	Position of which quantile should be used as a left (lower) boundary when calculating interquantile width.  Default value \code{NULL} results in using the start coordinate of the cluster.  Used only when \code{returnInterquantileWidth = TRUE}, otherwise ignored.
}
  \item{qUp}{
	Position of which quantile should be used as a right (upper) boundary when calculating interquantile width.  Default value \code{NULL} results in using the end coordinate of the cluster.  Used only when \code{returnInterquantileWidth = TRUE}, otherwise ignored.
}
}
\value{
	Returns a \code{data.frame} with genomic coordinates, position of dominant TSS, total CAGE signal and additional information for all TCs from specified CAGE dataset (sample).  If \code{returnInterquantileWidth = TRUE}, interquantile width for each TC is also calculated using specified quantile positions and returned in the data frame.
}
\author{Vanja Haberle}
\seealso{
\code{\link{clusterCTSS}}
}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

TC <- tagClusters(object = exampleCAGEset, sample = "sample2", returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)
head(TC)
}
