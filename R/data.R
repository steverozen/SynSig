#' Real exposure (signature attributions) from SignatureAnalyzer and SigProfiler
#'
#' @format Numerical matrix with rows indicating signatures and columns indicating
#' (tumor) samples.
#'
#' @note Prefix \code{sa} indicates SignatureAnalyzers, \code{sp} indicates
#' SigProfiler; \code{all} indicates all samples, \code{no.hyper} means
#' hypermutated tumours as defined for SignatureAnalyzer have been removed.
#'
#' @name RealExposures
NULL


#' @rdname RealExposures
#' @source \url{https://dx.doi.org/10.7303/syn11761237.4}
"sa.all.real.exposures"

#' @rdname RealExposures
"sp.all.real.exposures"

#' @rdname RealExposures
#' @source \url{https://dx.doi.org/10.7303/syn11761198.4}
"sa.no.hyper.real.exposures"

#' @rdname RealExposures
"sp.no.hyper.real.exposures"

#' Reference mutational signatures from PCAWG7
#'
#' @format Numerical matrix with rows indicating mutation
#' types and columns indicating signatures.
#'
#' @name MutationalSignatures

NULL

#' @rdname MutationalSignatures
"sa.96.sigs"

#' @rdname MutationalSignatures
"sa.COMPOSITE.sigs"

#' @rdname MutationalSignatures
"sp.sigs"

