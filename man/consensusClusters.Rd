\name{consensusClusters}
\alias{consensusClusters}
\alias{consensusClusters,CAGEset-method}
\title{Extracting consensus clusters from CAGEset object}
\description{
Extracts the information on consensus clusters from a CAGEset object.
}
\usage{
consensusClusters(object, sample = NULL, returnInterquantileWidth = FALSE, qLow = NULL, qUp = NULL)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
  \item{sample}{
	Optional. Label of the CAGE dataset (experiment, sample) for which to extract sample-specific information on consensus clusters.
}
  \item{returnInterquantileWidth}{
	Should the interquantile width of consensus clusters in specified sample be returned.  Used only when \code{sample} argument is specified, otherwise ignored.
}
  \item{qLow}{
	Position of which quantile should be used as a left (lower) boundary when calculating interquantile width.  Used only when \code{sample} argument is specified and \code{returnInterquantileWidth = TRUE}, otherwise ignored.
}
  \item{qUp}{
	Position of which quantile should be used as a right (upper) boundary when calculating interquantile width.  Used only when \code{sample} argument is specified and \code{returnInterquantileWidth = TRUE}, otherwise ignored.
}

}
\value{
	Returns a \code{data.frame} with information on consensus clusters, including genomic coordinates.  When \code{sample} argument is NOT specified, total CAGE signal across all CAGE datasets (samples) is returned in the \code{tpm} column.  When \code{sample} argument is specified, the \code{tpm} column contains CAGE signal of consensus clusters in that specific sample.  When \code{returnInterquantileWidth = TRUE}, additional sample-specific information is returned, including position of the dominant TSS, and interquantile width of the consensus clusters in the specified sample.
}
\author{Vanja Haberle}

\seealso{
\code{\link{tagClusters}}
}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

clusters.general <- consensusClusters(exampleCAGEset)
head(clusters.general)

clusters.sample <- consensusClusters(exampleCAGEset, sample = "sample2")
head(clusters.sample)
}
