\name{consensusClustersTpm}
\alias{consensusClustersTpm}
\alias{consensusClustersTpm,CAGEset-method}
\title{Extracting consensus clusters tpm matrix from CAGEset object}
\description{
Extracts the matrix with normalized CAGE tag values for consensus clusters across all samples from a CAGEset object.
}
\usage{
consensusClustersTpm(object)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
}
\value{
	Returns a \code{matrix} with normalized CAGE tag values across all samples.
}
\author{Vanja Haberle}

\seealso{
\code{\link{consensusClusters}}
}
\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

clusters.tpm <- consensusClustersTpm(exampleCAGEset)
head(clusters.tpm)
}
