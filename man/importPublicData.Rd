\name{importPublicData}
\alias{importPublicData}
\alias{importPublicData,list,character,character-method}
\title{Importing publicly available CAGE data from R data packages}
\description{
	Imports CAGE data from a dataset available in one of the two CAGE data packages: \code{FANTOM3and4CAGE} or \code{ENCODEprojectCAGE}, into a \code{\link{CAGEset}} object.  After the CAGEset object has been created the data can be further manipulated and visualized using other functions available in the \code{CAGEr} package.
}
\usage{
importPublicData(dataset, group, sample)
}
\arguments{
  \item{dataset}{
	One of the datasets available in \code{FANTOM3and4CAGE} or \code{ENCODEprojectCAGE} data packages.  It is a named \code{list} of \code{data.frames}, which contain genomic coordinates of CAGE detected TSSs and their frequency of usage.
}
  \item{group}{
	Character string specifying one or more groups from which the samples should be selected.  \code{group} must contain one or more names of the elements in the list provided as \code{dataset}.  Either only one group has to be specified (if all samples belong to the same group) or one group per sample (if samples belong to different groups).  In the latter case, the number of elements in \code{group} must match the number of elements in \code{sample}.  Check the corresponding data package for available groups in the specified dataset.
}
  \item{sample}{
	Character string specifying one or more CAGE samples.  Check the corresponding data package for available samples within each group and their labels.
}
}
\details{
	CAGE data produced by two big consortia (FANTOM and ENCODE) is available in R through two data packages: \code{FANTOM3and4CAGE} or \code{ENCODEprojectCAGE}, respectively.  \code{dataset} argument accepts any of the datasets from these packages.  \code{group} and \code{sample} arguments have to be the labels of the available groups/samples in the provided dataset.  If all specified samples belong to the same group, only this one group can be provided, otherwise, for each sample a corresponding group has to be specified, \emph{i.e.} the number of elements in \code{group} must match the numer of elements in \code{sample}.
}
\value{
	A \code{\link{CAGEset}} object is returned. Slots \code{librarySizes}, \code{CTSScoordinates} and \code{tagCountMatrix} are occupied by the data imported from the provided dataset for specified CAGE samples.
}
\author{Vanja Haberle}

\seealso{
\code{\link{getCTSS}}
}
\examples{
library(FANTOM3and4CAGE)
data(FANTOMmouseSamples)
FANTOMmouseSamples

data(FANTOMtissueCAGEmouse)

library(BSgenome.Mmusculus.UCSC.mm9)

exampleCAGEset <- importPublicData(dataset = FANTOMtissueCAGEmouse,
group = c("liver", "liver", "brain"), 
sample = c("cloned_mouse", "control_mouse", "brain"))

exampleCAGEset
}
