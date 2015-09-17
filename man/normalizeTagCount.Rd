\name{normalizeTagCount}
\alias{normalizeTagCount}
\alias{normalizeTagCount,CAGEset-method}
\title{Normalizing raw CAGE tag count}
\description{
Normalizes raw CAGE tag count per CTSS in all experiments to a same referent distribution.  A simple tag per million normalization or normalization to a referent power-law distribution (Balwierz \emph{et al}., Genome Biology 2009) can be specified.
}
\usage{
normalizeTagCount(object, method = "powerLaw", fitInRange = c(10, 1000),
                  alpha = 1.25, T = 10^6)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
  \item{method}{
	Method to be used for normalization.  Can be either \code{"simpleTpm"} to convert tag counts to tags per million or \code{"powerLaw"} to normalize to a referent power-law distribution, or \code{"none"} to keep using the raw tag counts in downstream analyses.
}
  \item{fitInRange}{
	An integer vector with two values specifying a range of tag count values to be used for fitting a power-law distribution to reverse cumulatives.  Used only when \code{method = "powerLaw"}, otherwise ignored.  See Details.
}
  \item{alpha}{
	\code{-1 * alpha} will be the slope of the referent power-law distribution in the log-log representation.  Used only when \code{method = "powerLaw"}, otherwise ignored.  See Details.
}
  \item{T}{
	Total number of CAGE tags in the referent power-law distribution.  Setting \code{T = 10^6} results in normalized values that correspond to tags per million in the referent distribution.  Used only when \code{method = "powerLaw"}, otherwise ignored.  See Details.
}
}
\details{
	It has been shown that many CAGE datasets follow a power-law distribution (Balwierz \emph{et al}., Genome Biology 2009).  Plotting the number of CAGE tags (X-axis) against the number of TSSs that are supported by >= of that number of tags (Y-axis) results in a distribution that can be approximated by a power-law.  On a log-log scale this theoretical referent distribution can be described by a monotonically decreasing linear function \code{y = -1 * alpha * x + beta}, which is fully determined by the slope \code{alpha} and total number of tags \code{T} (which together with \code{alpha} determines the value of \code{beta}).  Thus, by specifying parameters \code{alpha} and \code{T} a desired referent power-law distribution can be selected.  However, real CAGE datasets deviate from the power-law in the areas of very low and very high number of tags, so it is advisable to discard these areas before fitting a power-law distribution.  \code{fitInRange} parameter allows to specify a range of values (lower and upper limit of the number of CAGE tags) that will be used to fit a power-law.  Plotting reverse cumulatives using \code{\link{plotReverseCumulatives}} function can help in choosing the best range of values.  After fitting a power-law distribution to each CAGE dataset individually, all datasets are normalized to a referent distribution specified by \code{alpha} and \code{T}. When \code{T = 10^6}, normalized values are expressed as tags per million (tpm).
}
\value{
The slot \code{normalizedTpmMatrix} of the provided \code{\link{CAGEset}} object will be occupied by normalized CAGE signal values per CTSS across all experiments, or with the raw tag counts (in case \code{method = "none"}).
}
\references{
Balwierz \emph{et al}. (2009) Methods for analyzing deep sequencing expression data: constructing the human and mouse promoterome with deepCAGE data, \emph{Genome Biology} \bold{10}(7):R79. 
}
\author{Vanja Haberle}

\seealso{
\code{\link{plotReverseCumulatives}}
}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

normalizeTagCount(exampleCAGEset, method = "powerLaw")
}
