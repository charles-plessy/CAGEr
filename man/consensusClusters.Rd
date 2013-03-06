\name{consensusClusters}
\alias{consensusClusters}
\alias{consensusClusters,CAGEset-method}
\title{Extracting consensus clusters from CAGEset object}
\description{
Extracts the consensus clusters from a CAGEset object.
}
\usage{
consensusClusters(object)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
}
\value{
	Returns a \code{data.frame} with information on consensus clusters, including genomic coordinates and total CAGE signal across all CAGE datasets.
}
\author{Vanja Haberle}

\seealso{
\code{\link{tagClusters}}
}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

clusters <- consensusClusters(exampleCAGEset)
head(clusters)
}
