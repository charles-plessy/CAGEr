\name{setColors}
\alias{setColors}
\alias{setColors,CAGEset-method}
\title{Setting colors for samples}
\description{
	Assigns one color to each sample in the CAGEset object.  These colors are used in various plots and exported tracks to consistently represent corresponding samples.
}
\usage{
setColors(object, colors = NULL)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
	}
  \item{colors}{
	A character vector of valid color names.  For a complete list of valid color names see colors().  Alternatively, it can be a character vector of colors specified in hexadecimal format (\emph{e.g.} "#FF0000" for red).  Number of provided colors must match the number of samples in the CAGEset object.  Provided colors are assigned to samples according to their ordering in the CAGEset object, \emph{i.e} in the order they are returned by \code{\link{sampleLabels}} function.
}
}
\value{
Assigns one color to each sample in the CAGEset object by setting them as a name attribute of the \code{sampleLabels} slot.
}
\author{Vanja Haberle}

\seealso{
\code{\link{sampleLabels}}
}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

setColors(exampleCAGEset, colors = c("darkred", "navy", "forestgreen"))
sampleLabels(exampleCAGEset)
}
