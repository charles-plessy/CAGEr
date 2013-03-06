\name{CAGEr-package}
\alias{CAGEr-package}
\alias{CAGEr}
\docType{package}
\title{Analysis of CAGE (Cap Analysis of Gene Expression) sequencing data for precise mapping of transcription start sites and promoterome mining}
\description{
CAGEr package performs identification of transcription start sites and frequency of their usage from input CAGE sequencing data, normalization of raw CAGE tag count, clustering of TSSs into tag clusters (TC) and their aggregation across multiple CAGE experiments to construct the promoterome.  It manipulates multiple CAGE experiments at once, performs expression profiling across experiments both at level of individual TSSs and clusters of TSSs, exports several different types of track files for visualization in the UCSC Genome Browser, performs analysis of promoter width and detects differential usage of TSSs (promoter shifting) between samples. Multicore option for parallel processing is supported on Unix-like platforms.
}
\details{
\tabular{ll}{
Package: \tab CAGEr\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-12-18\cr
License: \tab GPL-3\cr
Depends: \tab R (>= 2.15.0), methods, BSgenome\cr
}
}
\author{
Vanja Haberle

Maintainer: Vanja Haberle <vanja.haberle@bio.uib.no>
}
\keyword{ package }
