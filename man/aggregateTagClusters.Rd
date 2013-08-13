\name{aggregateTagClusters}
\alias{aggregateTagClusters}
\alias{aggregateTagClusters,CAGEset-method}
\title{Aggregating tag clusters across multiple CAGE datasets}
\description{
	Aggregates tag clusters (TCs) across all CAGE dataset within the CAGEset object to create a referent set of consensus clusters.
}
\usage{
aggregateTagClusters(object, tpmThreshold = 5, 
					 excludeSignalBelowThreshold = TRUE, 
					 qLow = NULL, qUp = NULL, maxDist = 100)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
  \item{tpmThreshold}{
	Only tag clusteres with normalized signal \code{>= tpmThreshold} will be used to construct consensus clusters.
}
  \item{excludeSignalBelowThreshold}{
	When \code{TRUE} only tag clusteres with normalized signal \code{>= tpmThreshold} will contribute to the total CAGE signal of a consensus cluster, \emph{i.e.} only the TCs that are used to construct consensus cluster. When set to \code{FALSE} all TCs that overlap consensus cluster will contribute to the total signal (regardless whether they pass the threshold or not), however only the TCs above the threshold will be used to define consensus cluster boundaries. Thus, it that case the TCs above the threshold are first used to construct consensus clusters and define their boundaries, but then CAGE signal from all TCs that fall within those boundaries is used to calculate total signal of a particular consensus cluster.
}
  \item{qLow}{
	Position of which "lower" quantile should be used as 5' boundary of the tag cluster.  If \code{qLow = NULL} start position of the TC is used.  See Details.
}
  \item{qUp}{
	Position of which "upper" quantile should be used as 3' boundary of the tag cluster.  If \code{qUp = NULL} end position of the TC is used. \code{qUp} has to be \code{>= qLow}.  See Details.
}
  \item{maxDist}{
	Maximal length of the gap (in base-pairs) between two tag clusters for them to be part of the same consensus clusters.  See Details.
}
}
\details{
	Tag clusters (TCs) returned by \code{\link{clusterCTSS}} function are constructed for every CAGE dataset within CAGEset object separatelly, based on the CAGE signal in that sample. Thus, TCs from two CAGE datasets can differ both in their number, genomic coordinates, position of dominant TSS and overall signal. To be able to compare all samples at the level of clusters of TSSs, TCs from all CAGE datasets are aggregated into a single set of consensus clusters. First, TCs with signal \code{>= tpmThreshold} from all CAGE datasets are selected, and their 5' and 3' boundaries are determined based on provided \code{qLow} and \code{qUp} parameters. If \code{qLow = NULL} and \code{qUp = NULL} the start and end coordinates, \emph{i.e.} the full span of the TC is used, otherwise the positions of \code{qLow} and \code{qUp} quantiles are used as 5' and 3' boundary, respectively.  Finally, the defined set of TCs from all CAGE datasets is reduced to a non-overlapping set of consensus clusters by merging overlapping TCs and TCs \code{<= maxDist} base-pairs apart.  Consensus clusters represent a referent set of promoters that can be further used for expression profiling or detecting "shifting" (differentially used) promoters between different CAGE samples.
}
\value{
	The slots \code{consensusClusters}, \code{tagClustersInConsensusClusters} and \code{consensusClustersTpmMatrix} of the provided \code{\link{CAGEset}} object will be occupied by the genomic coordinates of consensus clusters, information on containing TCs and the total CAGE signal across all CAGE datasets, respectively.
}
\author{Vanja Haberle}
\seealso{
\code{\link{clusterCTSS}}
}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

aggregateTagClusters(object = exampleCAGEset, tpmThreshold = 50,
excludeSignalBelowThreshold = FALSE, qLow = 0.1, qUp = 0.9, maxDist = 100)
}
