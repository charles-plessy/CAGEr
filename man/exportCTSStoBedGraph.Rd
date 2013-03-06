\name{exportCTSStoBedGraph}
\alias{exportCTSStoBedGraph}
\alias{exportCTSStoBedGraph,CAGEset-method}
\title{Creating bedGraph tracks of CAGE transcription starts sites}
\description{
Creates bedGraph file(s) with track(s) of CAGE signal supporting each TSS that can be visualised in the UCSC Genome Browser.
}
\usage{
exportCTSStoBedGraph(object, values, oneFile = TRUE)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
  \item{values}{
	Specifies which values will be exported to the bedGraph file. Can be either \code{"raw"} to export raw tag count values or \code{"normalized"} to export normalized values.
}
  \item{oneFile}{
	Logical, should all CAGE datasets be exported as individual tracks into the same bedGraph file (TRUE) or into separate bedGraph files (FALSE).
}
}
\value{
	Creates bedGraph file(s) in the working directory that can be directly visualised as custom tracks in the UCSC Genome Browser.  If \code{oneFile = TRUE} one bedGraph file containing multiple annotated tracks will be created, otherwise two files per CAGE dataset will be created, one for plus strand and one for minus strand CTSSs, and they will be named according to the labels of individual datasets. 
}
\author{Vanja Haberle}

\seealso{
\code{\link{normalizeTagCount}}
}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

exportCTSStoBedGraph(exampleCAGEset, values = "normalized", oneFile = TRUE)
}
