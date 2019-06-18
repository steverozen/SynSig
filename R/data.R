#' Real exposure (signature attributions) from SignatureAnalyzer and SigProfiler
#'
#' @format Numerical matrix with rows indicating signatures and columns indicating
#' (tumor) samples.
#'
#' @note Prefix \code{sa} indicates SignatureAnalyzers, \code{sp} indicates
#' SigProfiler; \code{all} indicates all samples, \code{no.hyper} means
#' that hypermutated tumors as defined for SignatureAnalyzer have been removed.
#'
#' @name RealExposures
NULL


#' @rdname RealExposures
#' @source \url{https://dx.doi.org/10.7303/syn11761237.4}
"sa.all.real.exposures"

#' @rdname RealExposures
#' @source \url{https://dx.doi.org/10.7303/syn11738669.5}
"sp.all.real.exposures"

#' @rdname RealExposures
#' @source \url{https://dx.doi.org/10.7303/syn11761198.4}
"sa.no.hyper.real.exposures"

#' @rdname RealExposures
#' @source \url{https://dx.doi.org/10.7303/syn11761237.4}
"sp.no.hyper.real.exposures"

#' Reference mutational signature profiles from PCAWG7.
#'
#' @details \code{sa.96.sigs}
#' provides SignatureAnalyzer mutational signature profiles collapsed from
#' COMPOSITE to 96-channel SNS signatures.
#'
#' @details \code{sa.COMPOSITE.sigs} provides
#' COMPOSITE mutational signature profiles extracted by
#' SignatureAnalyzer. \code{sa.COMPOSITE.sigs}
#' are an \code{rbind} of the contents of
#' \url{https://www.synapse.org/#!Synapse:syn11738311} (SBS 1536),
#' \url{https://www.synapse.org/#!Synapse:syn11738308} (DBS), and
#' \url{https://www.synapse.org/#!Synapse:syn11738309} (ID).
#'
#' @details \code{sa.DBS.sigs} provides the DBS signatures extracted by
#' SignatureAnalyzer, from
#' \url{https://www.synapse.org/#!Synapse:syn11738312}. These are not the
#' DBS signatures that are part of \code{sa.COMPOSITE.sigs};
#' these were extracted from the ID catalogs alone.
#'
#' @details \code{sa.ID.sigs} provides the ID signatures extracted by
#' SignatureAnalyzer, from
#' \url{https://www.synapse.org/#!Synapse:syn11738313}.These are not the
#' ID signatures that are part of \code{sa.COMPOSITE.sigs};
#' these were extracted from the ID catalogs alone.
#'
#' @details \code{sp.sigs} provides signatures extracted by SigProfiler.
#'
#' @format Numerical matrix with rows indicating mutation
#' types and columns indicating signatures.
#'
#' @name MutationalSignatures

NULL

#' @rdname MutationalSignatures
#' @source \url{https://www.synapse.org/#!Synapse:syn11738310}
"sa.96.sigs"

#' @rdname MutationalSignatures
#' @source \url{https://www.synapse.org/#!Synapse:syn11738311}
#' @source \url{https://www.synapse.org/#!Synapse:syn11738308}
#' @source \url{https://www.synapse.org/#!Synapse:syn11738309}
"sa.COMPOSITE.sigs"

#' @rdname MutationalSignatures
#' @source \url{https://www.synapse.org/#!Synapse:syn11738312}
"sa.DBS.sigs"

#' @rdname MutationalSignatures
#' @source \url{https://www.synapse.org/#!Synapse:syn11738313}
"sa.ID.sigs"

#' @rdname MutationalSignatures
#' @source \url{https://www.synapse.org/#!Synapse:syn11738319}
"sp.sigs"


# Quiets concerns of R CMD check about no visible binding for global variable
if(getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("sa.all.real.exposures",
      "sp.all.real.exposures",
      "sa.COMPOSITE.sigs",
      "sa.96.sigs",
      "sp.sigs",
      "TEMPORARY",
      "BayesNMF.L1.KL.fixed_W.Z",
      "BayesNMF.L1.KL.fixed_W.Z.sample",
      "OutDir.dir"
      ))
}
