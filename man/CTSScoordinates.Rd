\name{CTSScoordinates}
\alias{CTSScoordinates}
\alias{CTSScoordinates,CAGEset-method}
\title{Extracting genomic coordinates of TSSs from CAGEset object}
\description{
Extracts the genomic coordinates of all detected TSSs from a CAGEset object.
}
\usage{
CTSScoordinates(object)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
}
\value{
	Returns a \code{data.frame} with genomic coordinates of all TSSs. \code{pos} column contains 1-based coordinate of the TSS.
}
\author{Vanja Haberle}

\seealso{
\code{\link{CTSStagCount}}
\cr
\code{\link{CTSSnormalizedTpm}}
}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

CTSS <- CTSScoordinates(exampleCAGEset)
head(CTSS)
}
