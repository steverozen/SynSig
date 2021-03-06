% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MapSPToSASignatureNamesInExposure.R
\name{MapSPToSASignatureNamesInExposure}
\alias{MapSPToSASignatureNamesInExposure}
\title{With the signatures represented in a matrix of exposures, find the nearest
SignatureAnalyzer exposure.}
\usage{
MapSPToSASignatureNamesInExposure(sp.exposures,
  sa.sig.names.to.consider = colnames(sa.96.sigs))
}
\arguments{
\item{sp.exposures}{The exposures}

\item{sa.sig.names.to.consider}{A subset of the colnames
of \code{\link{sa.96.sigs}}}
}
\value{
A list with

\enumerate{
\item \code{exp2} Copy of \code{sp.exposures} with the
rownames(signature names) updated according to the
match.

\item \code{sp.to.sa.sig.match}

\item \code{sa.to.sp.sig.match} Best matches in the opposite direction
}
}
\description{
With the signatures represented in a matrix of exposures, find the nearest
SignatureAnalyzer exposure.
}
\details{
IMPORTANT: uses the package global
variables \code{\link{sa.96.sigs}}
and \code{\link{sp.sigs}}.
}
