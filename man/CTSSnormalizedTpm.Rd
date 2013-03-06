\name{CTSSnormalizedTpm}
\alias{CTSSnormalizedTpm}
\alias{CTSSnormalizedTpm,CAGEset-method}
\title{Extracting normalized CAGE signal for TSSs from CAGEset object}
\description{
Extracts the normalized CAGE signal for all detected TSSs in all CAGE datasets from a CAGEset object.
}
\usage{
CTSSnormalizedTpm(object)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
}
\value{
	Returns a \code{data.frame} with normalized CAGE signal supporting each TSS (rows) in every CAGE dataset (columns).
}
\author{Vanja Haberle}

\seealso{
\code{\link{CTSScoordinates}}
\cr
\code{\link{CTSStagCount}}
}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

CAGEsignal <- CTSSnormalizedTpm(exampleCAGEset)
head(CAGEsignal)
}
