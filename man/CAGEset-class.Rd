\name{CAGEset-class}
\Rdversion{1.1}
\docType{class}
\alias{CAGEset-class}
\alias{CAGEset}

% Coercion:
\alias{coerce,data.frame,CAGEset-method}


\title{Class \code{"CAGEset"}}
\description{
This class is used to store one or more CAGE (Cap Analysis of Gene Expression) datasets in the form of TSSs derived from CAGE tags and frequency of their usage, and to store and extract all information generated during the workflow.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CAGEset", ...)}.
\cr
Objects of the class contain information on the genomic coordinates of TSSs derived from sequenced CAGE tags from multiple experiments, number of tags supporting each TSS in each experiment, normalized CAGE signal at each TSS, all information on specified parameters and results of all downstream analyses.  Object has to be created before reading in the data by specifying input files and their type, referent genome and labels for individual CAGE datasets, as described in the vignette.  Data is read by applying a function to a created object and all further slots are filled during the workflow by applying specific functions.
}
\section{Slots}{
  \describe{
    \item{\code{genomeName}:}{Object of class \code{"character"}: the name of the BSgenome package used as the referent genome }
    \item{\code{inputFiles}:}{Object of class \code{"character"}: the paths to input files }
    \item{\code{inputFilesType}:}{Object of class \code{"character"}: the type of input files (e.g. bam) }
    \item{\code{sampleLabels}:}{Object of class \code{"character"}: the labels of individual CAGE experiments }
    \item{\code{librarySizes}:}{Object of class \code{"integer"}: the total number of CAGE tags per experiment }
    \item{\code{CTSScoordinates}:}{Object of class \code{"data.frame"}: the genomic coordinates of CAGE transcription start sites (CTSSs) }
    \item{\code{tagCountMatrix}:}{Object of class \code{"data.frame"}: the number of CAGE tags supporting every CTSS in each experiment }
    \item{\code{normalizedTpmMatrix}:}{Object of class \code{"data.frame"}: the normalized CAGE signal supporting every CTSS in each experiment }
    \item{\code{CTSSexpressionClusteringMethod}:}{Object of class \code{"character"}: the method used for expression clustering of CTSSs}
    \item{\code{CTSSexpressionClasses}:}{Object of class \code{"character"}: the labels of expression classes of CTSSs returned by expression clustering }
    \item{\code{clusteringMethod}:}{Object of class \code{"character"}: the method used for clustering CTSSs into tag clusters (TC) }
    \item{\code{filteredCTSSidx}:}{Object of class \code{"logical"}: the index of CTSSs included in tag clusters }
    \item{\code{tagClusters}:}{Object of class \code{"list"}: the list of tag clusters per CAGE experiment }
    \item{\code{CTSScumulativesTagClusters}:}{Object of class \code{"list"}: the cumulative distribution of CAGE signal along TCs }
    \item{\code{tagClustersQuantileLow}:}{Object of class \code{"list"}: the positions of lower quantile(s) within TCs }
    \item{\code{tagClustersQuantileUp}:}{Object of class \code{"list"}: the positions of upper quantile(s) within TCs }
    \item{\code{tagClustersInConsensusClusters}:}{Object of class \code{"data.frame"}: the information on which TCs from which experiments are contained within each consensus cluster }
    \item{\code{consensusClusters}:}{Object of class \code{"data.frame"}: the genomic coordinates of consensus clusters created by aggregating TCs across experiments }
    \item{\code{consensusClustersTpmMatrix}:}{Object of class \code{"matrix"} the normalized CAGE signal for every consensus cluster in each experiment }
    \item{\code{consensusClustersExpressionClusteringMethod}:}{Object of class \code{"character"} the method used for expression clustering of consensus clusters }
    \item{\code{consensusClustersExpressionClasses}:}{Object of class \code{"character"} the labels of expression classes of consensus clusters returned by expression clustering }
    \item{\code{CTSScumulativesConsensusClusters}:}{Object of class \code{"list"} the cumulative distribution of CAGE signal along consensus clusters }
    \item{\code{consensusClustersQuantileLow}:}{Object of class \code{"list"} the positions of lower quantile(s) within consensus clusters }
    \item{\code{consensusClustersQuantileUp}:}{Object of class \code{"list"} the positions of upper quantile(s) within consensus clusters }
    \item{\code{shiftingGroupX}:}{Object of class \code{"character"} the label(s) of experiment(s) in the first shifting group }
    \item{\code{shiftingGroupY}:}{Object of class \code{"character"} the label(s) of experiment(s) in the second shifting group }
    \item{\code{consensusClustersShiftingScores}:}{Object of class \code{"data.frame"} the shifting scores and P-values/FDR for comparison of consensus clusters between two (groups of) experiments }
  }
}
\section{Methods}{
  \describe{
    \item{CTSSclusteringMethod}{\code{signature(object = "CAGEset")}: extracts the method used for clustering CTSSs into tag clusters (TC) }
    \item{CTSScoordinates}{\code{signature(object = "CAGEset")}: extracts the genomic coordinates of all CTSSs }
    \item{CTSSnormalizedTpm}{\code{signature(object = "CAGEset")}: extracts the normalized CAGE signal supporting every CTSS in each experiment }
    \item{CTSStagCount}{\code{signature(object = "CAGEset")}: extracts the number of CAGE tags supporting every CTSS in each experiment }
    \item{aggregateTagClusters}{\code{signature(object = "CAGEset")}: aggregates TCs across all experiments into consensus clusters }
    \item{clusterCTSS}{\code{signature(object = "CAGEset")}: clusters CTSSs into TCs per experiment }
    \item{consensusClusters}{\code{signature(object = "CAGEset")}: extracts the genomic coordinates and other information on consensus clusters }
    \item{consensusClustersTpm}{\code{signature(object = "CAGEset")}: extracts the matrix with tpm values for consensus clusters across all samples }
    \item{cumulativeCTSSdistribution}{\code{signature(object = "CAGEset")}: calculates the cumulative distribution of CAGE signal along TCs or consensus clusters }
    \item{exportCTSStoBedGraph}{\code{signature(object = "CAGEset")}: creates bedGraph files of CTSSs for visualization in the UCSC Genome Browser }
    \item{exportToBed}{\code{signature(object = "CAGEset")}: creates various types of BED files for visualization in the UCSC Genome Browser }
    \item{expressionClasses}{\code{signature(object = "CAGEset")}: extracts the labels of the expression classes of CTSSs or consensus clusters returned from expression profiling }
    \item{extractExpressionClass}{\code{signature(object = "CAGEset")}: extracts CTSSs or consensus clusters belonging to specified expression class }
    \item{genomeName}{\code{signature(object = "CAGEset")}: extracts the name of the BSgenome package used as the referent genome }
    \item{getCTSS}{\code{signature(object = "CAGEset")}: reads in specified input files and fills in information on detected CTSSs and their tag count }
    \item{getExpressionProfiles}{\code{signature(object = "CAGEset")}: performs expression clustering of CTSSs or consensus clusters across experiments }
    \item{getShiftingPromoters}{\code{signature(object = "CAGEset")}: extracts consensus clusters with shifting score and/or FDR above specified threshold }
    \item{inputFiles}{\code{signature(object = "CAGEset")}: extracts the paths of input CAGE data files }
    \item{inputFilesType}{\code{signature(object = "CAGEset")}: extracts the type of input CAGE data files }
    \item{librarySizes}{\code{signature(object = "CAGEset")}: extracts the library sizes of individual CAGE experiments within CAGEset object }
    \item{mergeSamples}{\code{signature(object = "CAGEset", mergeIndex = "numeric")}: merges specified experiments (samples) into one (e.g. replicas) }
    \item{normalizeTagCount}{\code{signature(object = "CAGEset")}: normalizes raw CAGE tag count }
    \item{plotCorrelation}{\code{signature(object = "CAGEset")}: plots pairwise scatter plots and calculates correlation between samples }
    \item{plotExpressionProfiles}{\code{signature(object = "CAGEset")}: creates file with beanplots of expression across experiments for CTSSs or consensus clusters belonging to different expression classes }
    \item{plotInterquantileWidth}{\code{signature(object = "CAGEset")}: creates file with histograms of interquantile width }
    \item{plotReverseCumulatives}{\code{signature(object = "CAGEset")}: creates file with reverse cumulative plots of CAGE tag count per CTSS }
    \item{quantilePositions}{\code{signature(object = "CAGEset")}: calculates the positions of specified quantiles within TCs or consensus clusters }
    \item{sampleLabels}{\code{signature(object = "CAGEset")}: extracts the labels of individual CAGE experiments within CAGEset object }
    \item{scoreShift}{\code{signature(object = "CAGEset", groupX = "character", groupY = "character")}: calculates the shifting score and tests the statistical significance of differential TSS usage for consensus clusters between two specified (groups of) samples }
    \item{setColors}{\code{signature(object = "CAGEset")}: assigns color to each sample to be used in visualisation }
    \item{show}{\code{signature(object = "CAGEset")}: displays CAGEset object in a user friendly way }
    \item{tagClusters}{\code{signature(object = "CAGEset")}: extracts the tag clusters for specified CAGE experiment }
	 }
}

\section{Coercion}{
 \describe{
    \item{}{
      \code{as(from, "CAGEset")}: Creates a CAGEset object from a data.frame object.
    }
  }
}

\author{Vanja Haberle}

\examples{
showClass("CAGEset")
}
\keyword{classes}
