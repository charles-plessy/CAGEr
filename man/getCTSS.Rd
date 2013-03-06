\name{getCTSS}
\alias{getCTSS}
\alias{getCTSS,CAGEset-method}
\title{Reading CAGE data from input file(s) and detecting TSSs}
\description{
	Reads input CAGE datasets into CAGEset object, constructs CAGE transcriptions start sites (CTSSs) and counts number of CAGE tags supporting every CTSS in each input experiment.  Preprocessing and quality filtering of input CAGE tags can be also performed before constructing TSSs.
}
\usage{
getCTSS(object, sequencingQualityThreshold = 10, 
        mappingQualityThreshold = 20, removeFirstG = TRUE)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
	}
  \item{sequencingQualityThreshold, mappingQualityThreshold}{
	Only CAGE tags with average sequencing quality \code{>= sequencingQualityThreshold} and mapping quality \code{>= mappingQualityThreshold} are kept. Used only if \code{inputFileType(object) == "bam"}, \emph{i.e} if input files are BAM files of aligned sequenced CAGE tags, otherwise ignored.  If there is no information on mapping quality in the BAM file (e.g. software used to align CAGE tags to the referent genome does not provide mapping quality) the \code{mappingQualityThreshold} parameter is ignored.
}
  \item{removeFirstG}{
	Logical, should the first nucleotide of the CAGE tag be removed in case it is a G and it does not map to the referent genome (\emph{i.e.} it is a mismatch).  See Details.
}
}
\details{
In the CAGE experimental protocol an additional G nucleotide is often attached to the 5' end of the tag by the template-free activity of the reverse transcriptase used to prepare cDNA (Harbers and Carninci, Nature Methods 2005).  In cases where there is a G at the 5' end of the CAGE tag that does not map to the corresponding genome sequence, it can confidently be considered spurious and should be removed from the tag to avoid misannotating actual TSS. Thus, setting \code{removeFirstG = TRUE} is highly recommended.
}
\value{
The slots \code{librarySizes}, \code{CTSScoordinates} and \code{tagCountMatrix} of the provided \code{\link{CAGEset}} object will be occupied by the information on CTSSs created from input CAGE files.
}
\references{
Harbers and Carninci (2005) Tag-based approaches for transcriptome research and genome annotation, \emph{Nature Methods} \bold{2}(7):495-502.
}
\author{Vanja Haberle}

\seealso{
\code{\link{CTSScoordinates}}
\cr
\code{\link{CTSStagCount}}
}
\examples{
library(BSgenome.Drerio.UCSC.danRer7)

pathsToInputFiles <- system.file("extdata", c("Zf.unfertilized.egg.chr17.ctss",
"Zf.30p.dome.chr17.ctss", "Zf.prim6.rep1.chr17.ctss"), package="CAGEr")
labels <- paste("sample", seq(1,3,1), sep = "")

myCAGEset <- new("CAGEset", genomeName = "BSgenome.Drerio.UCSC.danRer7",
inputFiles = pathsToInputFiles, inputFilesType = "ctss", sampleLabels = labels)

getCTSS(myCAGEset)
}
