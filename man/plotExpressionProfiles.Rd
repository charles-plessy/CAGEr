\name{plotExpressionProfiles}
\alias{plotExpressionProfiles}
\alias{plotExpressionProfiles,CAGEset-method}
\title{
	Plotting expression profiles derived from CAGE data
}
\description{
	Creates PDF file with beanplots representing distribution of normalized expression across CAGE experiments for individual expression classes.  Different expression classes are shown in different colors and are labeled according to the labels returned by expression clustering.
}
\usage{
plotExpressionProfiles(object, what)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
  \item{what}{
	Which level of expression clustering should be used.  Can be either \code{"CTSS"} to plot expression profiles of individual CTSSs or \code{"consensusClusters"} to plot expression profiles of consensus clusters.	
}
}
\details{
The created file contains beanplots representing distribution of normalized expression across CAGE experiments for individual expression classes shown in separate boxes.  Each labeled box represents one expression class and contains a set of beanplots - one per CAGE experiment. Individual CAGE experiments are shown on X-axis and scaled normalized expression on Y-axis.  Individual beanplots show distribution of normalized expression values of elements belonging to specific expression class in particular CAGE experiment, and the entire box represents single expression profile. Different expression classes (boxes) are plotted in different colors and are labeled with labels returned by expression clustering.
}
\value{
	Creates PDF file named "CTSS_expression_profiles.pdf" (in case \code{what = "CTSS"}) or "consensusClusters_expression_profiles.pdf" (in case \code{what = "consensusClusters"}) in the working directory.
}
\author{Vanja Haberle}
\seealso{
\code{\link{getExpressionProfiles}}
\cr
\code{\link{expressionClasses}}
\cr
\code{\link{extractExpressionClass}}
}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

plotExpressionProfiles(object = exampleCAGEset, what = "CTSS")
}
