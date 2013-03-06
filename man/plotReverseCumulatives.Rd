\name{plotReverseCumulatives}
\alias{plotReverseCumulatives}
\alias{plotReverseCumulatives,CAGEset-method}
\title{Plotting reverse cumulative number of CAGE tags per CTSS}
\description{
	Creates PDF file with plots of reverse cumulative number of CAGE tags per CTSS for all CAGE datasets present in the CAGEset object.  The plots should be used as help in choosing the parameters for power-law normalization: the range of values to fit the power-law and the slope of the referent power-law distribution (Balwierz \emph{et al}., Genome Biology 2009).
}
\usage{
plotReverseCumulatives(object, onePlot = FALSE)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
  \item{onePlot}{
	Logical, should all CAGE datasets be plotted in the same plot (TRUE) or in separate plots (FALSE) within the same PDF file. 
}
}
\value{
	Creates PDF file named "CTSS_reverse_cumulatives_all_samples.pdf" in the working directory.  The file contains plots of reverse cumulative number of CAGE tags per CTSS for each CAGE dataset within CAGEset object.
}
\references{
Balwierz \emph{et al}. (2009) Methods for analyzing deep sequencing expression data: constructing the human and mouse promoterome with deepCAGE data, \emph{Genome Biology} \bold{10}(7):R79. 
}
\author{Vanja Haberle}

\seealso{
\code{\link{normalizeTagCount}}
}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

plotReverseCumulatives(exampleCAGEset, onePlot = TRUE)
}
