\name{inputFiles}
\alias{inputFiles}
\alias{inputFiles,CAGEset-method}
\title{Extracting paths to input files from CAGEset object}
\description{
Extracts the paths to CAGE data input files from a CAGEset object.
}
\usage{
inputFiles(object)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
}
\value{
	Returns a character vector of paths to CAGE data input files.
}

\author{Vanja Haberle}

\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

inputFiles(exampleCAGEset)
}
