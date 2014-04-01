\name{exportToBed}
\alias{exportToBed}
\alias{exportToBed,CAGEset-method}
\title{Creating BED tracks of TSSs and clusters of TSSs}
\description{
	Creates BED file(s) with track(s) of individual CTSSs, tag clusters or consensus clusters. CTSSs and consensus clusters can be optionally colored in the color of their expression class. Tag clusters and consensus clusters can be displayed in a gene-like representation with a line showing full span on the cluster, filled block showing interquantile range and a thick box denoting position of the dominant (most frequently used) TSS.
}
\usage{
exportToBed(object, what, qLow = NULL, qUp = NULL, 
            colorByExpressionProfile = FALSE, oneFile = TRUE)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
  \item{what}{
	Which elements should be exported to BED track. Can be \code{"CTSS"} to export individual CTSSs, \code{"tagClusters"} to export tag clusters or \code{"consensusClusters"} to export consensus clusters.  See Details.
}
  \item{qLow}{
	Position of which "lower" quantile should be used as 5' boundary of the filled block in gene-like representation of the cluster.  Default value \code{NULL} uses start position of the cluster.  Ignored when \code{what = "CTSS"}.  See Details.
}
  \item{qUp}{
	Position of which "upper" quantile should be used as 3' boundary of the filled block in gene-like representation of the cluster.  Default value \code{NULL} uses end position of the cluster.  Ignored when \code{what = "CTSS"}.  See Details.
}
  \item{colorByExpressionProfile}{
	Logical, should blocks be colored in the color of their corresponding expression class. Ignored when \code{what = "tagClusters"}.
}
  \item{oneFile}{
	Logical, should all CAGE datasets be exported as individual tracks into the same BED file (TRUE) or into separate BED files (FALSE).  Ignored when \code{what = "CTSS"}, which by default produces only one track.
}
}
\details{
	This functions creates various representations of CTSSs, tag clusters and consensus clusters in the BED format, which can be directly visualised in the UCSC Genome Browser. 
	\cr
	\cr
	When \code{what = "CTSS"}, one BED file with single track of 1bp blocks representing all detected CTSSs (in all CAGE samples) is created.  CTSSs can be colored according to their expression class (provided the expression profiling of CTSSs was done by calling \code{\link{getExpressionProfiles}} function).  Colors of expression classes match the colors in which they are shown in the plot returned by the \code{\link{plotExpressionProfiles}} function.  For \code{colorByExpressionProfile = FALSE}, CTSSs included in the clusters are shown in black and CTSSs that were filtered out in gray. 
	\cr
	\cr
	When \code{what = "tagClusters"}, one track per CAGE dataset is created, which can be exported to a single BED file (by setting \code{oneFile = TRUE}) or separate BED files (by setting \code{oneFile = FALSE}).  In case \code{qLow = NULL} and \code{qUp = NULL}, TCs are represented as simple blocks showing only the full span of TC from  the start to the end.  Setting \code{qLow} and/or \code{qUp} parameters to a value of the desired quantile creates a gene-like representation with a line showing full span of the TC, filled block showing specified interquantile range and a thick 1bp block denoting position of the dominant (most frequently used) TSS.  All TCs in one track (one CAGE dataset) are shown in the same color.
	\cr
	\cr
	When \code{what = "consensusClusters"} consensus clusters are exported to BED file.  Since there is only one set of consensus clusters common to all CAGE datasets, only one track is created in case of a simple representation.  This means that when \code{qLow = NULL} and \code{qUp = NULL} one track with blocks showing the full span of consensus cluster from the start to the end is created. However, the distribution of the CAGE signal within consensus cluster can be different in different CAGE samples, resulting in different positions of quantiles and dominant TSS.  Thus, when \code{qLow} and/or \code{qUp} parameters are set to a value of the desired quantile, a separate track with a gene-like representation is created for every CAGE dataset.  These tracks can be exported to a single BED file (by setting \code{oneFile = TRUE}) or separate BED files (by setting \code{oneFile = FALSE}).  The gene-like representation is analogous to the one described above for the TCs.  In all cases consensus clusters can be colored according to their expression class (provided the expression profiling of consensus clusters was done by calling \code{\link{getExpressionProfiles}} function).  Colors of expression classes match the colors in which they are shown in the plot returned by the \code{\link{plotExpressionProfiles}} function.  For \code{colorByExpressionProfile = FALSE} all consensus clusters are shown in black.
	
}
\value{
	Creates BED file(s) in the working directory that can be directly visualised as custom tracks in the UCSC Genome Browser.
}
\author{Vanja Haberle}

\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

### exporting CTSSs colored by expression class
exportToBed(object = exampleCAGEset, what = "CTSS",
colorByExpressionProfile = TRUE)

### exporting tag clusters in gene-like representation
exportToBed(object = exampleCAGEset, what = "tagClusters",
qLow = 0.1, qUp = 0.9, oneFile = TRUE)
}
