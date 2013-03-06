\name{extractExpressionClass}
\alias{extractExpressionClass}
\alias{extractExpressionClass,CAGEset-method}
\title{Extracting elements of the specified expression class}
\description{
	Extracts CTSSs or consensus clusters belonging to a specified expression class.
}
\usage{
extractExpressionClass(object, what, which = "all")
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
  \item{what}{
	Which level of expression clustering should be used. Can be either \code{"CTSS"} to extract expression class of individual CTSSs or \code{"consensusClusters"} to extract expression class of consensus clusters.
}
  \item{which}{
	Which expression class should be extracted. It has to be one of the valid expression class labels (as returned by \code{\link{expressionClasses}} function), or \code{"all"} to extract members of all expression classes.
}
}
\value{
	Returns a \code{data.frame} of CTSSs (when \code{what = "CTSS"}) or consensus clusters (when \code{what = "consensusClusters"}) belonging to a specified expression class, with genomic coordinates and additional information.  Last column contains the label of the corresponding expression class.
}
\author{Vanja Haberle}

\seealso{
\code{\link{getExpressionProfiles}}
\cr
\code{\link{plotExpressionProfiles}}
\cr
\code{\link{expressionClasses}}
}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

CTSSexprClasses <- extractExpressionClass(exampleCAGEset, what = "CTSS",
which = "all")
head(CTSSexprClasses)
}
