\name{librarySizes}
\alias{librarySizes}
\alias{librarySizes,CAGEset-method}
\title{Extracting library sizes from CAGEset object}
\description{
Extracts the library sizes (total number of CAGE tags) for all CAGE datasets from a CAGEset object.
}
\usage{
librarySizes(object)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
}
\value{
	Returns an integer vector of total number of CAGE tags (library size) for all CAGE datasets in the CAGEset object.
}
\author{Vanja Haberle}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

librarySizes(exampleCAGEset)
}
