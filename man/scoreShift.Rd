\name{scoreShift}
\alias{scoreShift}
\alias{scoreShift,CAGEset,character,character-method}
\title{Calculating promoter shifting score}
\description{
Calculates shifting score for all consensus clusters (promoters) between two specified (groups of) CAGE datasets.  Shifting score is a measure of differential usage of TSSs within consensus cluster between two samples, which indicates the degree of physical separation of TSSs used in these samples within given consensus cluster.
}
\usage{
scoreShift(object, groupX, groupY, useMulticore = F, nrCores = NULL)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
  \item{groupX, groupY}{
	Character vector of the one or more CAGE dataset labels in the first (\code{groupX}) and in the second group (\code{groupY}).  Shifting score for each consensus cluster will be calculated by comparing CAGE signal in the samples from \code{groupX} against the signal in the samples from \code{groupY}.  If there is more than one CAGE dataset in the group, the datasets within that group will be merged together before comparison with the other group.  See Details.
}
  \item{useMulticore}{
	Logical, should multicore be used.  \code{useMulticore = TRUE} is supported only on Unix-like platforms.
}
  \item{nrCores}{
	Number of cores to use when \code{useMulticore = TRUE}.  Default value \code{NULL} uses all detected cores.
}
}
\details{
TSSs within one consensus cluster (promoter) can be used differently in different samples (cell types, tissues, developmental stages), with respect to their position and frequency of usage detected by CAGE.  This function calculates shifting scores of all consensus clusters between two specified (groups of) CAGE samples to detect promoters that are used differently in these two samples.  Shifting score is a measure of differential TSS usage defined as:
\cr
\cr
\code{score = max(F1 - F2) / max(F1)}
\cr
\cr
where F1 is a cumulative sum of CAGE signal along consensus cluster in the group of samples with lower total signal in that consensus cluster, and F2 in the opposite group.  Since cumulative sum can be calculated in both forward (5' -> 3') and reverse (3' -> 5') direction, shifting score is calculated for both cases and the bigger value is selected as final shifting score.  Value of the shifting score is in the range \code{[-Inf, 1]}, where value of \code{1} means complete physical separation of TSSs used in the two samples for given consensus cluster.  In general, any non-negative value of the shifting score can be interpreted as the proportion of transcription initiation in the sample with lower expression that is happening "outside" (either upstream or downstream) of the region used for transcription initiation in the other sample.  Negative values indicate no physical separation, \emph{i.e.} the region used for transcription initiation in the sample with lower expression is completely contained within the region used for transcription initiation in the other sample.
\cr
Since calculation of shifting score requires cumulative sums along all consensus clusters, they have to be calculated beforehand by calling \code{\link{cumulativeCTSSdistribution}} function.
}
\value{
The slots \code{shiftingGroupX}, \code{shiftingGroupY} and \code{consensusClustersShiftingScores} of the provided \code{\link{CAGEset}} object will be occupied by the information on the groups of CAGE datasets that have been compared and shifting scores of all consensus clusters.  Consensus clusters (promoters) with score above specified threshold can be extracted by calling \code{\link{getShiftingPromoters}} function.
}
\author{Vanja Haberle}

\seealso{
\code{\link{cumulativeCTSSdistribution}}
\cr
\code{\link{getShiftingPromoters}}
}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

scoreShift(object = exampleCAGEset, groupX = c("sample1", "sample2"),
groupY = "sample3")
}
