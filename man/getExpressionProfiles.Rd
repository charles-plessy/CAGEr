\name{getExpressionProfiles}
\alias{getExpressionProfiles}
\alias{getExpressionProfiles,CAGEset-method}
\title{
	CAGE data based expression clustering
}
\description{
	Performs clustering of CAGE derived expression across multiple experiments, both at level of individual TSSs or entire clusters of TSSs.
}
\usage{
getExpressionProfiles(object, what, tpmThreshold = 5, 
                      nrPassThreshold = 1, method = "som", 
                      xDim = 5, yDim = 5)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
  \item{what}{
	At which level should the expression clustering be done.  Can be either \code{"CTSS"} to perform clustering of individual CTSSs or \code{"consensusClusters"} to perform clustering of consensus clusters.  See Details.
}
  \item{tpmThreshold, nrPassThreshold}{
	Only CTSSs or consensus clusters (depending on \code{what} parameter) with normalized CAGE signal \code{>= tpmThreshold} in \code{>= nrPassThreshold} experiments will be included in expression clustering.
}
  \item{method}{
	Method to be used for expression clustering.  Can be either \code{"som"} to use the self-organizing map (SOM) algorithm (Toronen \emph{et al}., FEBS Letters 1999) implemented in the the \code{som} function from \code{som} package, or \code{"kmeans"} to use the K-means algorithm implemented in the \code{kmeans} function from \code{stats} package.
}
  \item{xDim, yDim}{
	When \code{method = "kmeans"}, \code{xDim} specifies number of clusters that will be returned by K-means algorithm and \code{yDim} is ignored.  When \code{method = "som"}, \code{xDim} specifies the the first and \code{yDim} the second dimension of the self-organizing map, which results in total \code{xDim * yDim} clusters returned by SOM.
}
}
\details{
	Expression clustering can be done at level of individual CTSSs, in which case the feature vector used as input for clustering algorithm contains log-transformed and scaled (divided by standard deviation) normalized CAGE signal at individual TSS across multiple experiments.  Only TSSs with normalized CAGE signal \code{>= tpmThreshold} in at least \code{nrPassThreshold} CAGE experiments are used for expression clustering.  However, CTSSs along the genome can be spatially clustered into tag clusters for each experiment separately using the \code{\link{clusterCTSS}} function, and then aggregated across experiments into consensus clusters using \code{\link{aggregateTagClusters}} function.  Once the consensus clusters have been created, expression clustering at the level of these wider genomic regions (representing entire promoters rather than individual TSSs) can be performed.  In that case the feature vector used as input for clustering algorithm contains normalized CAGE signal within entire consensus cluster across multiple experiments, and threshold values in \code{tpmThreshold} and \code{nrPassThreshold} are applied to entire consensus clusters.
}
\value{
	If \code{what = "CTSS"} the slots \code{CTSSexpressionClusteringMethod} and \code{CTSSexpressionClasses} will be occupied, and if \code{what = "consensusClusters"} the slots \code{consensusClustersExpressionClusteringMethod} and \code{consensusClustersExpressionClasses} of the provided \code{\link{CAGEset}} object will be occupied with the results of expression clustering.  Labels of expression classes (clusters) can be retrieved using \code{\link{expressionClasses}} function, and elements belonging to a specific expression class can be selected using \code{\link{extractExpressionClass}} function.
}
\references{
Toronen \emph{et al}. (1999) Analysis of gene expression data using self-organizing maps, \emph{FEBS Letters} \bold{451}:142-146.
}
\author{Vanja Haberle}

\seealso{
\code{\link{plotExpressionProfiles}}
\cr
\code{\link{expressionClasses}}
\cr
\code{\link{extractExpressionClass}}
}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

getExpressionProfiles(object = exampleCAGEset, what = "CTSS",
tpmThreshold = 50, nrPassThreshold = 1, method = "som", xDim = 3, yDim = 3)
}
