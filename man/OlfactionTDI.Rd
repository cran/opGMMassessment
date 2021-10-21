\name{OlfactionTDI}
\alias{OlfactionTDI}

\docType{data}
\title{Example data of olfactory subtest performance.}

\description{Data set containing the results of 3 odor subtests collected from 1,000 subjects and assigned to the clinical olfactory diagnostic categories of anosmia, hyposmia, or normosmia.}

\usage{data("OlfactionTDI")}

\details{
Size 1000 x 3 , stored in \code{OlfactionTDI$[Threshold, Discrimination, Identification]}
Original classes 3, stored in \code{OlfactionTDI$Class}
}

\examples{
data(OlfactionTDI)
str(OlfactionTDI)
}
