\name{inputFilesType}
\alias{inputFilesType}
\alias{inputFilesType,CAGEset-method}
\title{Extracting type of input files from CAGEset object}
\description{
Extracts the information on the type of CAGE data input files from a CAGEset object.
}
\usage{
inputFilesType(object)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
}
\value{
	Returns the label of the file type of CAGE data input files, \emph{e.g.} \code{"bam"} or \code{"ctss"}.
}
\author{Vanja Haberle}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

inputFilesType(exampleCAGEset)
}
