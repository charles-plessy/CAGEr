\name{importPublicData}
\alias{importPublicData}
\alias{importPublicData,character,character,ANY,character-method}
\title{Importing publicly available CAGE data from various resources}
\description{
	Imports CAGE data from different sources into a \code{\link{CAGEset}} object.  After the \code{CAGEset} object has been created the data can be further manipulated and visualized using other functions available in the \emph{CAGEr} package and integrated with other analyses in R.  Available resources include:
	\cr
	- FANTOM5 datasets (Forrest \emph{et al.} Nature 2014) for numerous human and mouse samples (primary cells, cell lines and tissues), which are fetched directly from FANTOM5 online resource.
	\cr
	- FANTOM3 and 4 datasets (Carninci \emph{et al.} Science 2005, Faulkner \emph{et al.} Nature Genetics 2009, Suzuki \emph{et al.} Nature Genetics 2009) from \emph{FANTOM3and4CAGE} data package available from Bioconductor
	\cr
	- ENCODE datasets (Djebali \emph{et al.} Nature 2012) for numerous human cell lines from \emph{ENCODEprojectCAGE} data package, which is available for download from \href{http://promshift.genereg.net/CAGEr/}{http://promshift.genereg.net/CAGEr/}.
	\cr
	- Zebrafish developmental timecourse datasets (Nepal \emph{et al.} Genome Research 2013) from \emph{ZebrafishDevelopmentalCAGE} data package, which is available for download from \href{http://promshift.genereg.net/CAGEr/}{http://promshift.genereg.net/CAGEr/}.
}
\usage{
importPublicData(source, dataset, group, sample)
}
\arguments{
  \item{source}{
	Character vector specifying one of the available resources for CAGE data. Can be one of the following:
	\cr
	\code{"FANTOM5"}: for fetching and importing CAGE data for various human or mouse primary cells, cell lines and tissues from the online FANTOM5 resource (\href{http://fantom.gsc.riken.jp/5/data/}{http://fantom.gsc.riken.jp/5/data/}).  All data published in main FANTOM5 publication by Forrest \emph{et al.} is available.
	\cr
	\code{"FANTOM3and4"}: for importing CAGE data for various human or mouse tissues produced within FANTOM3 and FANTOM4 projects. Requires data package \emph{FANTOM3and4CAGE} to be installed. This data package is available from Bioconductor.
	\cr
	\code{"ENCODE"}: for importing CAGE data for human cell lines from ENCODE project published by Djebali \emph{et al.}. Requires data package \emph{ENCODEprojectCAGE} to be installed. This data package is available for download from \href{http://promshift.genereg.net/CAGEr/}{http://promshift.genereg.net/CAGEr/}.
	\cr
	\code{"ZebrafishDevelopment"}: for importing CAGE data from developmental timecourse of zebrafish (\emph{Danio rerio}) published by Nepal \emph{et al.}. Requires data package \emph{ZebrafishDevelopmentalCAGE} to be installed. This data package is available for download from \href{http://promshift.genereg.net/CAGEr/}{http://promshift.genereg.net/CAGEr/}.
	\cr
	See Details for further explanation of individual resources.
}
  \item{dataset}{
	Character vector specifying one or more of the datasets available in the selected resource. For \code{FANTOM5} it can be either \code{"human"} or \code{"mouse"}, and only one of them can be specified at a time. For other resources please refer to the vignette of the corresponding data package for the list of available datasets. Multiple datasets mapped to the same genome can be specified to combine selected samples from each.
}
  \item{group}{
	Character string specifying one or more groups within specified dataset(s), from which the samples should be selected.  \code{group} argument is used only when importing TSSs from data packages and ignored when \code{source="FANTOM5"}.  For available groups in each dataset please refer to the vignette of the corresponding data package. Either only one group has to be specified (if all selected samples belong to the same group) or one group per sample (if samples belong to different groups).  In the latter case, the number of elements in \code{group} must match the number of elements in \code{sample}.
}
  \item{sample}{
	Character string specifying one or more CAGE samples.  Check the corresponding data package for available samples within each group and their labels. For FANTOM5 resource, list of all human (~1000) and mouse (~) samples can be obtained in \emph{CAGEr} by loading \code{data(FANTOM5humanSamples)} and \code{data(FANTOM5mouseSamples)}, respectively.  Use the names from the \code{sample} column to specify which samples should be imported.
}
}
\details{
	CAGE data from different sources is available for importing directly into \code{CAGEset} object for further manipulation with \emph{CAGEr}.
	\cr
	FANTOM consortium provides single base-pair resolution TSS data for numerous human and mouse primary cells, cell lines and tissues produced within FANTOM5 project (Forrest \emph{et al.} Nature 2014).  These are directly fetched from their online resource at \href{http://fantom.gsc.riken.jp/5/data}{http://fantom.gsc.riken.jp/5/data} and imported into a \code{CAGEset} object.  To use this resource specify \code{source="FANTOM5"}.  The \code{dataset} argument can be either \code{"human"} or \code{"mouse"}, but not both at the same time.  The list of all human and mouse samples can be obtained by loading \code{data(FANTOM5humanSamples)} and \code{data(FANTOM5mouseSamples)}.  The \code{sample} column gives the names of individual samples that should be provided as \code{sample} argument.  See example below.
	\cr
	TSS data from previous FANTOM3 and FANTOM4 projects (Carninci \emph{et al.}, Faulkner \emph{et al.}, Suzuki \emph{et al.}) are also available through \emph{FANTOM3and4CAGE} data package.  This data package can be installed directly from Bioconductor.  To use this resource install and load \emph{FANTOM3and4CAGE} package and specify \code{source="FANTOM3and4"}. The \code{dataset} argument can be a name of any of the datasets available in this package.  Load \code{data(FANTOMhumanSamples)} or \code{data(FANTOMmouseSamples)} for the list of available datasets with group and sample labels for specific human or mouse samples.  These have to be provided as \code{dataset}, \code{group} and \code{sample} arguments to import selected samples.  If all samples belong to the same group, only this one group can be provided, otherwise, for each sample a corresponding group has to be specified, \emph{i.e.} the number of elements in \code{group} must match the numer of elements in \code{sample}.
	\cr
	ENCODE consortium produced CAGE data for numerous human cell lines (Djebali \emph{et al.} Nature 2012). We have used these data to derive single base-pair resolution TSSs and collected them into an R data package \emph{ENCODEprojectCAGE}.  This data package is available for download from \href{http://promshift.genereg.net/CAGEr/}{http://promshift.genereg.net/CAGEr/}.  To use this resource install and load \emph{ENCODEprojectCAGE} data package and specify \code{source="ENCODE"}.  The \code{dataset} argument can be a name of any of the datasets available in this package.  Load \code{data(ENCODEhumanCellLinesSamples)} for the list of available datasets with group and sample labels for specific samples.  These have to be provided as \code{dataset}, \code{group} and \code{sample} arguments to import selected samples.  Multiple datasets can be combined together, by specifying them as \code{dataset} argument.  If all samples belong to the same dataset and the same group, these dataset and group can be specified only once, otherwise, for each sample a corresponding dataset and group have to be specified, \emph{i.e.} the number of elements in \code{dataset} and \code{group} must match the numer of elements in \code{sample}.
	\cr
	Precise TSSs are also available for zebrafish (\emph{Danio Rerio}) from CAGE data published by Nepal \emph{et al.} for 12 developmental stages.  These have been collected into a data package \emph{ZebrafishDevelopmentalCAGE}, which is available for download from \href{http://promshift.genereg.net/CAGEr/}{http://promshift.genereg.net/CAGEr/}.  To use this resource install and load \emph{ZebrafishDevelopmentalCAGE} data package and specify \code{source="ZebrafishDevelopment"}.  Load \code{data(ZebrafishSamples)} for the list of available datasets and group and sample labels, which have to be specified to import these data.
}
\value{
	A \code{\link{CAGEset}} object is returned. Slots \code{librarySizes}, \code{CTSScoordinates} and \code{tagCountMatrix} are occupied by the single base-pair resolution TSS data imported from the selected resource.
}
\references{
	Carninci \emph{et al.} (2005) The Transcriptional Landscape of the Mammalian Genome, \emph{Science} \bold{309}(5740):1559-1563.
	\cr
	Djebali \emph{et al.} (2012) Landscape of transcription in human cells, \emph{Nature} \bold{488}(7414):101-108.
	\cr
	Faulkner \emph{et al.} (2009) The regulated retrotransposon transcriptome of mammalian cells, \emph{Nature Genetics} \bold{41}:563-571.
	\cr
	Forrest \emph{et al.} (2014) A promoter-level mammalian expression atlas, \emph{Nature} \bold{507}(7493):462-470.
	\cr
	Nepal \emph{et al.} (2013) Dynamic regulation of the transcription initiation landscape at single nucleotide resolution during vertebrate embryogenesis, \emph{Genome Research} \bold{23}(11):1938-1950.
	\cr
	Suzuki \emph{et al.} (2009) The transcriptional network that controls growth arrest and differentiation in a human myeloid leukemia cell line, \emph{Nature Genetics} \bold{41}:553-562.
}
\author{Vanja Haberle}

\seealso{
\code{\link{getCTSS}}
}
\examples{
### importing FANTOM5 data

# list of FANTOM5 human tissue samples
data(FANTOM5humanSamples)
head(subset(FANTOM5humanSamples, type == "tissue"))

# import selected samples
exampleCAGEset <- importPublicData(source="FANTOM5", dataset = "human", sample = c("adipose_tissue__adult__pool1", "adrenal_gland__adult__pool1", "aorta__adult__pool1"))

exampleCAGEset


### importing FANTOM3/4 data from a data package
library(FANTOM3and4CAGE)

# list of mouse datasets available in this package
data(FANTOMmouseSamples)
unique(FANTOMmouseSamples$dataset)
head(subset(FANTOMmouseSamples, dataset == "FANTOMtissueCAGEmouse"))
head(subset(FANTOMmouseSamples, dataset == "FANTOMtimecourseCAGEmouse"))

# import selected samples from two different mouse datasets
exampleCAGEset <- importPublicData(source="FANTOM3and4", dataset = c("FANTOMtissueCAGEmouse", "FANTOMtimecourseCAGEmouse"), group = c("brain", "adipogenic_induction"), sample = c("CCL-131_Neuro-2a_treatment_for_6hr_with_MPP+", "DFAT-D1_preadipocytes_2days"))

exampleCAGEset
}
