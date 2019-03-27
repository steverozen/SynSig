#' SynSig: Create catalogs of synthetic mutational spectra and assess the
#' performance of mutational-signature analysis programs on these.
#'
#' @section Overview:
#'
#' The main focus is generating synthetic catalogs of mutational
#' spectra (mutations in tumors) based on known mutational signature
#' profiles and attributions (assignment of exposures to tumors) in
#' the PCAWG7 data. We call this kind of synthetic
#' data broadly "reality-based" synthetic
#' data.  The package also has a set of functions that
#' generate random mutational signature profiles and then create
#' synthetic catalogs based on these random signature profiles. We
#' call this kind of synthetic data "random" synthetic data, while
#' pointing out that much depends on the distributions from which
#' the random signature profiles and attributions are generated.
#'
#' Typical workflow for generating catalogs of "reality-based" synthetic
#' mutational spectra is as follows.
#'
#' \preformatted{
#'
#' Input (based on SignatureAnalyzer or SigProfiler analysis of PCAWG tumors)
#'   A, matrix of attributions (signatures x samples)
#'   S, mutational signature profiles (mutation type x signature)
#'
#' P <- GetSynSigParamsFromExposures(A, ...)
#'
#' synthetic.exposures <- GenerateSyntheticExposures(P, ...)
#'
#' synthetic.spectra <- CreateAndWRiteCatalog(S, synthetic.exposures, ...)
#'
#' T <- Signatures extracted by SignatureAnalzer or SigProfiler on synthetic.spectra
#'
#' SummarizeResults(T, S, synthetic.exposures, ...)
#'
#' }
#'
#' @section Creating Synthetic Mutational Catalogs:
#'
#' These functions create synthetic mutational catalogs based
#' on parameters derived from signature profiles
#' and attributions (exposures).
#'
#' @section Summarize results (of signature extraction):
#'
#' Relevant functions are: \enumerate{
#'
#' \item \code{\link{SummarizeSigProfiler}}
#' \item \code{\link{SignatureAnalyzerSummarizeTopLevel}}
#' \item \code{\link{SignatureAnalyzerSummarizeSBS1SBS5}}
#'
#' }
#'
#' @section Comparing two sets of mutational signatures:
#'
#' Functions for comparing mutational signatures and
#' sets of mutational signatures. Often we will be interested
#' in comparing signature profiles extracted from synthetic data to the
#' ground-truth signature profiles.
#'
#' \code{\link{Match1Sig}},
#' \code{\link{MatchSigs1Direction}},
#' \code{\link{MatchSigs2Directions}},
#' \code{\link{MatchSigsAndRelabel}}
#'
#' @docType package
#' @name SynSig

NULL
