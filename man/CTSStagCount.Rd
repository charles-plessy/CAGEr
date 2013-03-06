\name{CTSStagCount}
\alias{CTSStagCount}
\alias{CTSStagCount,CAGEset-method}
\title{Extracting CAGE tag count for TSSs from CAGEset object}
\description{
Extracts the tag count for all detected TSSs in all CAGE datasets from a CAGEset object.
}
\usage{
CTSStagCount(object)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
}
\value{
	Returns a \code{data.frame} with number of CAGE tags supporting each TSS (rows) in every CAGE dataset (columns).
}
\author{Vanja Haberle}

\seealso{
\code{\link{CTSScoordinates}}
\cr
\code{\link{CTSSnormalizedTpm}}
}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

tagCount <- CTSStagCount(exampleCAGEset)
head(tagCount)
}
