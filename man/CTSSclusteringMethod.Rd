\name{CTSSclusteringMethod}
\alias{CTSSclusteringMethod}
\alias{CTSSclusteringMethod,CAGEset-method}
\title{Extracting CTSS clustering method from CAGEset object}
\description{
Extracts the label of the method used for CTSS clustering into tag clusters from a CAGEset object.
}
\usage{
CTSSclusteringMethod(object)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
}
\value{
	Returns a label of the method used for CTSS clustering.
}
\author{Vanja Haberle}

\seealso{
\code{\link{clusterCTSS}}
}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

CTSSclusteringMethod(exampleCAGEset)
}
