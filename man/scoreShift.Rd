\name{scoreShift}
\alias{scoreShift}
\alias{scoreShift,CAGEset,character,character-method}
\title{Calculating promoter shifting score}
\description{
Calculates shifting score for all consensus clusters (promoters) between two specified (groups of) CAGE datasets.  Shifting score is a measure of differential usage of TSSs within consensus cluster between two samples, which indicates the degree of physical separation of TSSs used in these samples within given consensus cluster. In addition to shifting score, a statistical significance (P-value and FDR) of differential TSS usage is calculated for each consensus cluster using Kolmogorov-Smirnov test.
}
\usage{
scoreShift(object, groupX, groupY, testKS = TRUE, useTpmKS = TRUE, useMulticore = F, nrCores = NULL)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
  \item{groupX, groupY}{
	Character vector of the one or more CAGE dataset labels in the first (\code{groupX}) and in the second group (\code{groupY}).  Shifting score for each consensus cluster will be calculated by comparing CAGE signal in the samples from \code{groupX} against the signal in the samples from \code{groupY}.  If there is more than one CAGE dataset in the group, the datasets within that group will be merged together before comparison with the other group.  See Details.
}
  \item{testKS}{
	Logical, should Kolomogorov-Smirnov test for statistical significance of differential TSS usage be performed, and P-values and FDR returned.  See Details.
}
  \item{useTpmKS}{
	Logical, should normalized (tpm) values (\code{TRUE}) or raw tag counts (\code{FALSE}) be used to derive sample sizes for Kolomogorov-Smirnov test.  Used only when \code{testKS = TRUE}, otherwise ignored.  See Details.
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
In addition to shifting score which indicates only physical separation (upstream or downstream shift of TSSs), a more general assessment of differential TSS usage can be obtained by performing a two-sample Kolmogorov-Smirnov test on cumulative sums of CAGE signal along the consensus cluster.  In that case, cumulative sums in both samples are scaled to range [0,1] and are considered to be empirical cumulative distribution functions (ECDF) reflecting sampling of TSS positions during transcription initiation.  Kolmogorov-Smirnov test is performed to assess whether the two underlying probability distributions differ.  To obtain P-value (\emph{i.e.} the level at which the null-hypothesis can be rejected), sample sizes that generated the ECDFs are required, in addition to actual K-S statistics calculated from ECDFs.  These are derived either from raw tag counts, \emph{i.e.} exact number of times each TSS in the cluster was sampled during sequencing (when \code{useTpmKS = FALSE}), or from normalized tpm values (when \code{useTpmKS = TRUE}). P-values obtained from K-S tests are further adjusted for multiple testing using Benjamini & Hochberg (BH) method and for each P-value a corresponding false-discovery rate (FDR) is also reported.
\cr
Since calculation of shifting scores and Kolmogorov-Smirnov test require cumulative sums along consensus clusters, they have to be calculated beforehand by calling \code{\link{cumulativeCTSSdistribution}} function.
}
\value{
The slots \code{shiftingGroupX}, \code{shiftingGroupY} and \code{consensusClustersShiftingScores} of the provided \code{\link{CAGEset}} object will be occupied by the information on the groups of CAGE datasets that have been compared and shifting scores of all consensus clusters.  Consensus clusters (promoters) with shifting score and/or FDR above specified threshold can be extracted by calling \code{\link{getShiftingPromoters}} function.
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
groupY = "sample3", testKS = TRUE, useTpmKS = FALSE)
}
