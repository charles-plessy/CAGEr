\name{plotInterquantileWidth}
\alias{plotInterquantileWidth}
\alias{plotInterquantileWidth,CAGEset-method}
\title{Plotting distribution of interquantile width}
\description{
	Creates PDF file with histograms showing distribution of interquantile width of tag clusters or consensus clusters in each CAGE dataset.
}
\usage{
plotInterquantileWidth(object, clusters, tpmThreshold = 5,
                       qLow = 0.1, qUp = 0.9, xlim = c(0,150), ...)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
  \item{clusters}{
	Which clusters should be used.  Can be either \code{"tagClusters"} to plot distribution of interquantile width of tag clusters or \code{"consensusClusters"} to plot distribution of interquantile width of consensus clusters.
}
  \item{tpmThreshold}{
	Only clusters with normalized signal \code{>= tpmThreshold} will be included in the histogram.
}
  \item{qLow}{
	Position of which "lower" quantile should be used as 5' boundary.  See Details.
}
  \item{qUp}{
	Position of which "upper" quantile should be used as 3' boundary.  See Details.
}
  \item{xlim}{
	The x axis limits of the plot.
}
  \item{...}{
	Additional arguments passed to \code{plot()} function, such as \code{ylim}, \emph{etc.}.
}
}
\details{
Interquantile width is a width (in base-pairs) of the central part of the genomic region (bounded by the positions of specified \code{qLow} and \code{qUp} quantiles) that contains \code{>= (qUp - qLow) * 100\%} of the CAGE signal.  Positions of specified quantiles within each cluster have to be calculated beforehand by calling \code{\link{quantilePositions}} function.  Interquantile width is a more robust measure of the promoter width than the total span of the region, because it takes into account the magnitude of the expression in the region.
}
\value{
Creates PDF file named "tagClusters_interquantile_width_all_samples.pdf" or "consensusClusters_interquantile_width_all_samples.pdf" in the working directory (depending on the value of \code{cluster} parameter).  The file contains histograms showing distribution of interquantile width in every CAGE experiment.
}
\author{Vanja Haberle}

\seealso{
\code{\link{quantilePositions}}
}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

plotInterquantileWidth(object = exampleCAGEset, clusters = "tagClusters",
tpmThreshold = 50, qLow = 0.1, qUp = 0.9)
}
