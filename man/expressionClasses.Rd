\name{expressionClasses}
\alias{expressionClasses}
\alias{expressionClasses,CAGEset-method}
\title{Extracting labels of expression classes}
\description{
Retrieves labels of expression classes of either individual CTSSs or consensus clusters from a CAGEset object.
}
\usage{
expressionClasses(object, what)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
  \item{what}{
	Which level of expression clustering should be used. Can be either \code{"CTSS"} to extract labels of expression classes of individual CTSSs or \code{"consensusClusters"} to extract labels of expression classes of consensus clusters.
}
}
\value{
	Returns character vector of labels of expression classes.  The number of labels matches the number of expression clusters returned by \code{\link{getExpressionProfiles}} function.
}
\author{Vanja Haberle}

\seealso{
	\code{\link{getExpressionProfiles}}
	\cr
	\code{\link{plotExpressionProfiles}}
	\cr
	\code{\link{extractExpressionClass}}
}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

exprClasses <- expressionClasses(exampleCAGEset, what = "CTSS")
exprClasses
}
