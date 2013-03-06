\name{clusterCTSS}
\docType{methods}
\alias{clusterCTSS}
\alias{clusterCTSS,CAGEset-method}
\title{Clustering CTSSs into tag clusters (TCs)}
\description{
Clusters individual CAGE transcription start sites (CTSSs) along the genome into tag clusters using specified "ab initio" method, or assigns them to predefined genomic regions.
}
\usage{
clusterCTSS(object, threshold = 1, nrPassThreshold = 1,
            thresholdIsTpm = TRUE, method = "distclu", maxDist = 20, 
            removeSingletons = FALSE, keepSingletonsAbove = Inf, 
            minStability = 1, maxLength = 500, 
            reduceToNonoverlapping = TRUE, customClusters = NULL, 
            useMulticore = FALSE, nrCores = NULL)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
  \item{threshold, nrPassThreshold}{
	Only CTSSs with signal \code{>= threshold} in \code{>= nrPassThreshold} experiments will be used for clustering and will contribute towards total signal of the cluster.
}
  \item{thresholdIsTpm}{
	Logical, is threshold raw tag count value (FALSE) or normalized signal (TRUE)
}
  \item{method}{
	Method to be used for clustering. Can be one of the \code{"distclu"}, \code{"paraclu"} or \code{"custom"}.  See Details.
}
  \item{maxDist}{
	Maximal distance between two neighbouring CTSSs for them to be part of the same cluster.  Used only when \code{method = "distclu"}, otherwise ignored.
}
  \item{removeSingletons}{
	Logical, should tag clusters containing only one CTSS be removed.  Ignored when \code{method = "custom"}.
}
  \item{keepSingletonsAbove}{
	Controls which singleton tag clusters will be removed.  When \code{removeSingletons = TRUE}, only singletons with signal \code{< keepSingletonsAbove} will be removed.  Useful to prevent removing highly supported singleton tag clusters.  Default value \code{Inf} results in removing all singleton TCs when \code{removeSingletons = TRUE}.  Ignored when \code{removeSingletons = FALSE} or \code{method = "custom"}.
}
  \item{minStability}{
	Minimal stability of the cluster, where stability is defined as ratio between maximal and minimal density value for which this cluster is maximal scoring.  For definition of stability refer to Frith \emph{et al}., Genome Research, 2007.  Clusters with stability \code{< minStability} will be discarded.  Used only when \code{method = "paraclu"}, otherwise ignored. 
}
  \item{maxLength}{
	Maximal length of cluster in base-pairs.  Clusters with length \code{> maxLength} will be discarded.  Ignored when \code{method = "custom"}.
}
  \item{reduceToNonoverlapping}{
	Logical, should smaller clusters contained within bigger cluster be removed to make a final set of tag clusters non-overlapping. Used only when \code{method = "paraclu"}.  See Details.
}
  \item{customClusters}{
	Genomic coordinates of predefined regions to be used to segment the CTSSs.  It has to be a \code{data.frame} with following columns: \code{chr} (chromosome name), \code{start} (0-based start coordinate), \code{end} (end coordinate), \code{strand} (either \code{"+"}, or \code{"-"}).  Used only when \code{method = "custom"}.
}
  \item{useMulticore}{
	Logical, should multicore be used.  \code{useMulticore = TRUE} is supported only on Unix-like platforms.
}
  \item{nrCores}{
	Number of cores to use when \code{useMulticore = TRUE}.  Default value \code{NULL} uses all detected cores.
}
}
\details{
Two "ab initio" methods for clustering TSSs along the genome are supported: \code{"distclu"} and \code{"paraclu"}.  \code{"distclu"} is an implementation of simple distance-based clustering of data attached to sequences, where two neighbouring TSSs are joined together if they are closer than some specified distance.  \code{"paraclu"} is an implementation of Paraclu algorithm for parametric clustering of data attached to sequences developed by M. Frith (Frith \emph{et al.}, Genome Research, 2007, \href{http://www.cbrc.jp/paraclu/}{http://www.cbrc.jp/paraclu/}).  Since Paraclu finds clusters within clusters (unlike distclu), additional parameters (\code{removeSingletons}, \code{keepSingletonsAbove}, \code{minStability}, \code{maxLength} and \code{reduceToNonoverlapping}) can be specified to simplify the output by discarding too small (singletons) or too big clusters, and to reduce the clusters to a final set of non-overlapping clusters.  Clustering is done for every CAGE dataset within CAGEset object separatelly, resulting in a different set of tag clusters for every CAGE dataset. TCs from different datasets can further be aggregated into a single referent set of consensus clusters by calling \code{\link{aggregateTagClusters}} function.
}
\value{
The slots \code{clusteringMethod}, \code{filteredCTSSidx} and \code{tagClusters} of the provided \code{\link{CAGEset}} object will be occupied by the information on method used for clustering, CTSSs included in the clusters and list of tag clusters per CAGE experiment, respectively.  To retrieve tag clusters for individual CAGE dataset use \code{\link{tagClusters}} function.
}
\references{
	Frith \emph{et al.} (2007) A code for transcription initiation in mammalian genomes, \emph{Genome Research} \bold{18}(1):1-12, (\href{http://www.cbrc.jp/paraclu/}{http://www.cbrc.jp/paraclu/}).
}
\author{Vanja Haberle}

\seealso{
\code{\link{tagClusters}}
\cr
\code{\link{aggregateTagClusters}}
}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

clusterCTSS(object = exampleCAGEset, threshold = 50, thresholdIsTpm = TRUE,
nrPassThreshold = 1, method = "distclu", maxDist = 20, 
removeSingletons = TRUE, keepSingletonsAbove = 100)
}
