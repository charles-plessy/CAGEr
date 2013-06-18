\name{plotCorrelation}
\alias{plotCorrelation}
\alias{plotCorrelation,CAGEset-method}
\title{
	Plotting pairwise scatter plots of CAGE tag counts per TSS and calculating correlation between samples 
}
\description{
	Creates PNG file with scatter plots of CAGE tag counts per TSS for all pairs of selected samples and calculates the correlation between them.
}
\usage{
plotCorrelation(object, samples = "all", method = "pearson",
				tagCountThreshold = 1, applyThresholdBoth = FALSE)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
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
}
\value{
	Creates PNG file named "Pairwise_tag_count_correlation.png" in the working directory.  Returns a \code{matrix} of pairwise correlations between selected samples.
}
\author{Vanja Haberle}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

plotCorrelation(object = exampleCAGEset)
}
