% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Create.syn.30.random.R
\name{CreateOneSetOfRandomCatalogs}
\alias{CreateOneSetOfRandomCatalogs}
\title{Create a pair of "random" synthetic catalogs, one for 96-channel
features and one for COMPOSITE features, for one set
of signatures.}
\usage{
CreateOneSetOfRandomCatalogs(num.syn.tumors, total.num.sigs, mut.mean,
  mut.sd, num.sigs.mean, num.sigs.sd, sig.name.prefix, sample.name.prefix,
  composite.dir.name, x96.dir.name, COMPOSITE.features,
  overwrite = FALSE)
}
\arguments{
\item{num.syn.tumors}{Total number of synthetic tumors to create.}

\item{total.num.sigs}{Total number of signatures in the universe.}

\item{mut.mean}{Mean of the log10 of the
number of mutations due to each signature.}

\item{mut.sd}{Standard deviation of the log10 of
the number of mutations due to each signature.}

\item{num.sigs.mean}{Mean number of signatures contributing to each tumor.}

\item{num.sigs.sd}{Standard deviation the number number of signatures
contribution to each tumor.}

\item{sig.name.prefix}{String to put in front of an integer (as
string) to form an identifier for a synthetic signature.}

\item{sample.name.prefix}{String to put in front of an integer (as
string) to form an identifier for a synthetic sample (tumor).}

\item{composite.dir.name}{string indicating the name of the COMPOSITE
subdirectory; probably one of \code{"sa.sa.COMPOSITE"} or
\code{"sp.sa.COMPOSITE"}.}

\item{x96.dir.name}{A string indicating the name of the 96-channel
subdirectory; probably one of \code{"sa.sa.COMPOSITE"} or
\code{"sp.sa.COMPOSITE"}.}

\item{COMPOSITE.features}{Character vector containing
rownames for a COMPOSITE signature or catalog.}

\item{overwrite}{If \code{TRUE} overwrite existing directories / files.}
}
\description{
Create a pair of "random" synthetic catalogs, one for 96-channel
features and one for COMPOSITE features, for one set
of signatures.
}
\keyword{internal}
