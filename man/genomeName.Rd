\name{genomeName}
\alias{genomeName}
\alias{genomeName,CAGEset-method}
\title{Extracting genome name from CAGEset object}
\description{
Extracts the name of a referent genome from a CAGEset object.
}
\usage{
genomeName(object)
}
\arguments{
  \item{object}{
	A \code{\link{CAGEset}} object
}
}
\value{
	Returns a name of a BSgenome package used as a referent genome.
}
\author{Vanja Haberle}

\examples{
load(system.file("data", "exampleCAGEset.RData", package="CAGEr"))

genomeName(exampleCAGEset)
}
