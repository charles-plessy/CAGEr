\name{plotCorrelation}
\alias{plotCorrelation}
\alias{plotCorrelation,CAGEset-method}
\title{
	Plotting pairwise scatter plots of CAGE signal and calculating correlation between samples 
}
\description{
	Creates PNG file with scatter plots of CAGE signal for all pairs of selected samples and calculates the correlation between them.
}
\usage{
plotCorrelation(object, what = "CTSS", values = "raw",
				samples = "all", method = "pearson",
				tagCountThreshold = 1, applyThresholdBoth = FALSE,
				plotSize = 800)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
  \item{what}{
	Which level should be used for plotting and calculating correlation.  Can be either \code{"CTSS"} to use individual TSSs or \code{"consensusClusters"} to use consensus clusters, \emph{i.e.} entire promoters.
}
  \item{values}{
	Specifies which values will be used for plotting and calculating correlation.  Can be either \code{"raw"} to use raw tag count per TSS or \code{"normalized"} to use normalized CAGE signal.  Used only when \code{what = "CTSS"}, otherwise ignored.
}
  \item{samples}{
	Character vector of sample labels for which the scatter plots should be plotted and correlation calculated.  Can be either \code{"all"} to plot and calculate pairwise correlations between all samples in a \code{\link{CAGEset}} object, or a subset of valid sample labels as returned by \code{\link{sampleLabels}} function.	
}
  \item{method}{
	A character string indicating which correlation coefficient should be computed.  Passed to \code{cor} function.  Can be one of \code{"pearson"}, \code{"spearman"}, or \code{"kendall"}.
}
  \item{tagCountThreshold, applyThresholdBoth}{
	Only TSSs with tag count \code{>= tagCountThreshold} in either one (\code{applyThresholdBoth = FALSE}) or both samples (\code{applyThresholdBoth = TRUE}) are plotted and used to calculate correlation.
}
  \item{plotSize}{
	Size of the individual comparison plot in pixels - the total size of the resulting png will be \code{length(samples) * plotSize} in both dimensions.
}
}
\value{
	Creates PNG file named "Pairwise_tag_count_correlation.png" in the working directory.  Returns a \code{matrix} of pairwise correlations between selected samples.
}
\author{Vanja Haberle}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

plotCorrelation(object = exampleCAGEset)
}
