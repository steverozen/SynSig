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
"sa.all.real.exposures"

#' @rdname RealExposures
"sp.all.real.exposures"

#' @rdname RealExposures
"sa.no.hyper.real.exposures"

#' @rdname RealExposures
"sp.no.hyper.real.exposures"

#' Signatures from SignatureAnalyzer
#'
#' @format Numerical matrix with rows indicating mutation
#' types and columns indicating signatures.
#'
#' @name SignatureAnalyzerSigs
NULL

#' @rdname SignatureAnalyzerSigs
"sa.96.sigs"

#' @rdname SignatureAnalyzerSigs
"sa.COMPOSITE.sigs"

