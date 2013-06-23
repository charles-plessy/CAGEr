\name{plotReverseCumulatives}
\alias{plotReverseCumulatives}
\alias{plotReverseCumulatives,CAGEset-method}
\title{Plotting reverse cumulative number of CAGE tags per CTSS}
\description{
	Creates PDF file with plots of reverse cumulative number of CAGE tags per CTSS for all CAGE datasets present in the CAGEset object.  The plots should be used as help in choosing the parameters for power-law normalization: the range of values to fit the power-law and the slope of the referent power-law distribution (Balwierz \emph{et al}., Genome Biology 2009).
}
\usage{
plotReverseCumulatives(object, values = "raw", 
				fitInRange = c(10,1000), onePlot = FALSE)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
  \item{values}{
	Specifies which values should be plotted.  Can be either \code{"raw"} to plot reverse cumulatives of raw CAGE tag counts or \code{"normalized"} to plot normalized tag count values.
}
  \item{fitInRange}{
	An integer vector with two values specifying a range of tag count values to be used for fitting a power-law distribution to reverse cumulatives.  Used only when \code{values = "raw"}, otherwise ignored.  See Details.
}
  \item{onePlot}{
	Logical, should all CAGE datasets be plotted in the same plot (TRUE) or in separate plots (FALSE) within the same PDF file. 
}
}
\details{
Number of CAGE tags (X-axis) is plotted against the number of TSSs that are supported by >= of that number of tags (Y-axis) on a log-log scale for each sample. In addition, a power-law distribution is fitted to each reverse cumulative using the values in the range specified by \code{fitInRange} parameter. The fitted distribution is defined by \code{y = -1 * alpha * x + beta} on the log-log scale, and the value of \code{alpha} for each sample is shown on the plot. In addition, a suggested referent power-law distribution to which all samples should be normalized is drawn on the plot and corresponding parameters (slope alpha and total number of tags T) are denoted on the plot.  Referent distribution is chosen so that its slope (alpha) is the median of slopes fitted to individual samples and its total number of tags (T) is the power of 10 nearest to the median number of tags of individual samples.  Resulting plots are helpful in deciding whether power-law normalization is appropriate for given samples and reported \code{alpha} values aid in choosing optimal \code{alpha} value for referent power-law distribution to which all samples will be normalized. For details about normalization see \code{\link{normalizeTagCount}} function.
}
\value{
	Creates PDF file named "CTSS_reverse_cumulatives_*_all_samples.pdf" in the working directory, where * denotes either "raw" or "normalized" depending on specified \code{values} parameter.  The file contains plots of reverse cumulative number of CAGE tags per CTSS for each CAGE dataset within CAGEset object.  Alpha values of fitted power-laws and suggested referent power-law distribution are reported on the plot in case \code{values = "raw"}.
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

plotReverseCumulatives(exampleCAGEset, values = "raw",
fitInRange = c(10,500), onePlot = TRUE)
}
