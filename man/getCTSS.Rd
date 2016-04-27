\name{getCTSS}
\alias{getCTSS}
\alias{getCTSS,CAGEset-method}
\title{Reading CAGE data from input file(s) and detecting TSSs}
\description{
	Reads input CAGE datasets into CAGEset object, constructs CAGE transcriptions start sites (CTSSs) and counts number of CAGE tags supporting every CTSS in each input experiment.  Preprocessing and quality filtering of input CAGE tags, as well as correction of CAGE-specific 'G' nucleotide addition bias can be also performed before constructing TSSs.
}
\usage{
getCTSS(object, sequencingQualityThreshold = 10, 
        mappingQualityThreshold = 20, removeFirstG = TRUE, 
        correctSystematicG = TRUE)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
	}
  \item{sequencingQualityThreshold, mappingQualityThreshold}{
	Only CAGE tags with average sequencing quality \code{>= sequencingQualityThreshold} and mapping quality \code{>= mappingQualityThreshold} are kept. Used only if \code{inputFileType(object) == "bam"} or \code{inputFileType(object) == "bamPairedEnd"}, \emph{i.e} when input files are BAM files of aligned sequenced CAGE tags, otherwise ignored.  If there are no sequencing quality values in the BAM file (\emph{e.g.} HeliScope single molecule sequencer does not return sequencing qualities) all reads will by default have this value set to -1.  Since the default value of \code{sequencingQualityThreshold} is 10, all the reads will consequently be discarded. To avoid this behaviour and keep all sequenced reads set \code{sequencingQualityThreshold} to -1 when processing data without sequencing qualities.  If there is no information on mapping quality in the BAM file (\emph{e.g.} software used to align CAGE tags to the referent genome does not provide mapping quality) the \code{mappingQualityThreshold} parameter is ignored. In case of paired-end sequencing BAM file (i.e. \code{inputFileType(object) == "bamPairedEnd"}) only the first mate of the properly paired reads (i.e. the five prime end read) will be read and subject to specified thresholds.
}
  \item{removeFirstG}{
	Logical, should the first nucleotide of the CAGE tag be removed in case it is a G and it does not map to the referent genome (\emph{i.e.} it is a mismatch).  Used only if \code{inputFileType(object) == "bam"} or \code{inputFileType(object) == "bamPairedEnd"}, \emph{i.e} when input files are BAM files of aligned sequenced CAGE tags, otherwise ignored.  See Details.
}
  \item{correctSystematicG}{
	Logical, should the systematic correction of the first G nucleotide be performed for the positions where there is a G in the CAGE tag and G in the genome.  This step is performed in addition to removing the first G of the CAGE tags when it is a mismatch, \emph{i.e.} this option can only be used when \code{removeFirstG = TRUE}, otherwise it is ignored.  The frequency of adding a G to CAGE tags is estimated from mismatch cases and used to systematically correct the G addition for positions with G in the genome.  Used only if \code{inputFileType(object) == "bam"} or \code{inputFileType(object) == "bamPairedEnd"}, \emph{i.e} when input files are BAM files of aligned sequenced CAGE tags, otherwise ignored.  See Details.
}
}
\details{
In the CAGE experimental protocol an additional G nucleotide is often attached to the 5' end of the tag by the template-free activity of the reverse transcriptase used to prepare cDNA (Harbers and Carninci, Nature Methods 2005).  In cases where there is a G at the 5' end of the CAGE tag that does not map to the corresponding genome sequence, it can confidently be considered spurious and should be removed from the tag to avoid misannotating actual TSS. Thus, setting \code{removeFirstG = TRUE} is highly recommended.
\cr
However, when there is a G both at the beginning of the CAGE tag and in the genome, it is not clear whether the original CAGE tag really starts at this position or the G nucleotide was added later in the experimental protocol.  To systematically correct CAGE tags mapping at such positions, a general frequency of adding a G to CAGE tags can be calculated from mismatch cases and applied to estimate the number of CAGE tags that have G added and should actually start at the next nucleotide/position.  The option \code{correctSystematicG} is an implementation of the correction algorithm described in Carninci \emph{et al.}, Nature Genetics 2006, Supplementary Information section 3-e.
}
\value{
The slots \code{librarySizes}, \code{CTSScoordinates} and \code{tagCountMatrix} of the provided \code{\link{CAGEset}} object will be occupied by the information on CTSSs created from input CAGE files.
}
\references{
Harbers and Carninci (2005) Tag-based approaches for transcriptome research and genome annotation, \emph{Nature Methods} \bold{2}(7):495-502.
\cr
Carninci \emph{et al.} (2006) Genome-wide analysis of mammalian promoter architecture and evolution, \emph{Nature Genetics} \bold{38}(7):626-635.
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
