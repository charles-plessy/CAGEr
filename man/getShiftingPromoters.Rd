\name{getShiftingPromoters}
\alias{getShiftingPromoters}
\alias{getShiftingPromoters,CAGEset-method}
\title{Selecting consensus clusters with shifting score above specified threshold}
\description{
Extracts consensus clusters with shifting score and/or FDR (adjusted P-value from Kolmogorov-Smirnov test) above specified threshold. Returns their genomic coordinates, total CAGE signal and the position of dominant TSS in the two compared groups of CAGE samples, along with the value of the shifting score, P-value and FDR.  Scores and P-values/FDR have to be calculated beforehand by calling \code{\link{scoreShift}} function.
}
\usage{
getShiftingPromoters(object, tpmThreshold = 0, scoreThreshold = -Inf, fdrThreshold = 1)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
  \item{tpmThreshold}{
	Consensus clusters with total CAGE signal \code{>= tpmThreshold} in each of the compared groups will be returned.
}
  \item{scoreThreshold}{
	Consensus clusters with shifting score \code{>= scoreThreshold} will be returned. The default value \code{-Inf} returns all consensus clusters (for which score could be calculated, \emph{i.e.} the ones that have at least one tag in each of the comapred samples).
}
  \item{fdrThreshold}{
	Consensus clusters with adjusted P-value (FDR) from Kolmogorov-Smirnov test \code{>= fdrThreshold} will be returned. The default value \code{1} returns all consensus clusters (for which K-S test could be performed, \emph{i.e.} the ones that have at least one tag in each of the compared samples).
}
}
\value{
	Returns a \code{data.frame} of shifting promoters with genomic coordinates and positions of dominant TSS and CAGE signal in the two compared (groups of) samples, along with shifting score and adjusted P-value (FDR).
}
\author{Vanja Haberle}

\seealso{
\code{\link{scoreShift}}
}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

shifting.promoters <- getShiftingPromoters(object = exampleCAGEset,
tpmThreshold = 100, scoreThreshold = 0.4, fdrThreshold = 0.01)
head(shifting.promoters)
}
