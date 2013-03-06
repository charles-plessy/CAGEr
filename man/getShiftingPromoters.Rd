\name{getShiftingPromoters}
\alias{getShiftingPromoters}
\alias{getShiftingPromoters,CAGEset-method}
\title{Selecting consensus clusters with shifting score above specified threshold}
\description{
Extracts consensus clusters with shifting score above specified threshold and returns their genomic coordinates, total CAGE signal and the position of dominant TSS in the two compared groups of CAGE samples, and the value of the shifting score.  Scores have to be calculated beforehand by calling \code{\link{scoreShift}} function.
}
\usage{
getShiftingPromoters(object, tpmThreshold, scoreThreshold)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
  \item{tpmThreshold}{
	Consensus clusters with total CAGE signal \code{>= tpmThreshold} in each of the compared groups will be returned.
}
  \item{scoreThreshold}{
	Consensus clusters with shifting score \code{>= scoreThreshold} will be returned.
}
}
\value{
	Returns a \code{data.frame} of shifting promoters with genomic coordinates and positions of dominant TSS and CAGE signal in the two compared (groups of) samples.
}
\author{Vanja Haberle}

\seealso{
\code{\link{scoreShift}}
}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

shifting.promoters <- getShiftingPromoters(object = exampleCAGEset,
tpmThreshold = 100, scoreThreshold = 0.6)
head(shifting.promoters)
}
