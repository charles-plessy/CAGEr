\name{FANTOM5humanSamples}
\docType{data}
\alias{FANTOM5humanSamples}
\title{List of FANTOM5 human CAGE samples}
\description{
This is a \code{data.frame} object that contains information about all human primary cell, cell line and tissue CAGE datasets (988 datasets) produced by FANTOM5 constortium, which are available for import with \code{CAGEr}.  The list containes all samples published in the main FANTOM5 publication (Forrest \emph{et al.} Nature 2014) as presented in the Supplementary Table 1.  The columns provide the following information:
\cr
\code{sample}: name/label of individual sample, which can be used as \code{sample} argument in the \code{\link{importPublicData}} function to retrieve specified sample.
\cr
\code{type}: type of the sample, either \code{"cell line"}, \code{"primary cell"} or \code{"tissue"}.
\cr
\code{description}: description of individual sample as provided by the FANTOM5 consortium.
\cr
\code{library_id}: unique ID of the CAGE library as provided by the FANTOM5 consortium.
\cr
\code{data_url}: URL to gzip-ed TSS file at online FANTOM5 data resource, which is used by \code{\link{importPublicData}} to fetch given sample.
}
\usage{data(FANTOM5humanSamples)}
\format{A \code{data.frame} object}
\keyword{datasets}