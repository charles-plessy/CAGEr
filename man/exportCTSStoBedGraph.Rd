\name{exportCTSStoBedGraph}
\alias{exportCTSStoBedGraph}
\alias{exportCTSStoBedGraph,CAGEset-method}
\title{Creating bedGraph/bigWig tracks of CAGE transcription starts sites}
\description{
Creates bedGraph or BigWig file(s) with track(s) of CAGE signal supporting each TSS that can be visualised in the UCSC Genome Browser.
}
\usage{
exportCTSStoBedGraph(object, values = "normalized", format = "BigWig", oneFile = TRUE)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
  \item{values}{
	Specifies which values will be exported to the bedGraph file. Can be either \code{"raw"} to export raw tag count values or \code{"normalized"} to export normalized values.
}
  \item{format}{
	The format of the output. 
}
  \item{oneFile}{
	Logical, should all CAGE datasets be exported as individual tracks into the same bedGraph file (TRUE) or into separate bedGraph files (FALSE). Used only when \code{format="bedGraph"}, otherwise ignored.
}
}
\value{
	Creates bedGraph or BigWig file(s) in the working directory that can be directly visualised as custom tracks in the UCSC Genome Browser.  If \code{format="bedGraph"} and \code{oneFile = TRUE} one bedGraph file containing multiple annotated tracks will be created, otherwise two files per CAGE dataset will be created, one for plus strand and one for minus strand CTSSs, and they will be named according to the labels of individual datasets.  All bedGraph files contain headers with track description and can be directly uploaded as custom tracks to the UCSC Genome Browser. 
\cr
When \code{format="bigWig"}, two binary BigWig files per CAGE dataset are created, one for plus strand and one for minus strand CTSSs. Since BigWig files cannot contain headers with track description, a separate file named "CTSS.normalized.all.samples.track.description.txt" is created, which contains track headers for all BigWig files. To use these headers for adding custom tracks to the UCSC Genome Browser, move the BigWig files to a web location and edit the bigDataUrl sections in the headers file to point to corresponding BigWig files.
}
\author{Vanja Haberle}

\seealso{
\code{\link{normalizeTagCount}}
}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

exportCTSStoBedGraph(exampleCAGEset, values = "normalized", format = "bedGraph", oneFile = TRUE)
}
