\name{sampleLabels}
\alias{sampleLabels}
\alias{sampleLabels,CAGEset-method}
\title{Extracting CAGE datasets labels from CAGEset object}
\description{
Extracts the labels of CAGE datasets (samples, experiments) from a CAGEset object.
}
\usage{
sampleLabels(object)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
}
\value{
	Returns a character vector of labels of all CAGE datasets present in the CAGEset object.
}
\author{Vanja Haberle}

\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

sampleLabels(exampleCAGEset)
}
