#' Generate parallel synthetic exposures from SA and SP attributions and signatures
#'
#' @param sa.real.exp Exposure matrix from SignatureAnalzyer.
#'
#' @param sp.realexp Exposure matrix from SignatureProfiler.
#'
#' @param ca.type The type the cancer, which is used in sample identifiers,
#'   which SigProfiler expects. For example if \code{ca.type}
#'   is \code{Kidney-RCC}, sample identifiers will be of the form
#'   \code{Kidney-RCC::S9999}.
#'
#' @param num.syn.tumors Number of synthetic tumors to generate.
#'
#' @param file.prefix A file prefix to add before the exposure and
#'  parameter file.
#'
#'  @return A list with the following elements:
#'  \enumerate{
#'
#'  \item \code{sa.parms} The paramters computed from \code{sa.real.exp}.
#'  This a matrix with a column for each signature
#'  and 3 rows: \enumerate{
#'    \item The proportion of tumors with
#'  a given signature (in \code{sa.real.exp})
#'  \item The mean of the log10 of the number of mutations
#'  \item The standard
#'  deviation of log10 of the number of mutations

#'  \item \code{sa.syn.exp} The synthetic exposures computed from
#'  \code{sa.parms}
#'  \item \code{sp.parms} The parameters computed from \code{sp.real.exp}.
#'  \item \code{sp.syn.exp} The synthetic exposures computed from
#'  \code{sp.parms}
#'
#'  }
#'
#'  @details Creates a bunch of files in location
#'  governed by \code{\link{OutDir}}. The main raison d'etre of the
#'  function is to ensure that some convetions regarding file
#'  naming are followed.
#'
#'  This function does \strong{not} create the synthetic
#'  mutational spectra catalogs.

SAAndSPSynDataOneCAType <-
  function(sa.real.exp,
           sp.real.exp,
           ca.type, # Used in sample identifiers
           num.syn.tumors,
           file.prefix) {
    ca.type <- paste0(ca.type, "::")
    samples.to.use <-
      grep(ca.type, colnames(sa.real.exp), fixed = TRUE)
    sa.exp <- sa.real.exp[ , samples.to.use]
    sp.exp <- sp.real.exp[ , samples.to.use]
    stopifnot(colnames(sa.exp) == colnames(sp.exp))
    sa.info <-
      GenerateSynFromReal(
        sa.exp, num.syn.tumors,
        file.prefix = paste0("sa.", file.prefix),
        sample.id.prefix = paste0("SA.Syn.", ca.type, "S"))
    sp.info <-
      GenerateSynFromReal(
        sp.exp, num.syn.tumors,
        file.prefix = paste0("sp.", file.prefix),
        sample.id.prefix = paste0("SP.Syn.", ca.type, "S"))
    return(list(
      sa.parms  = sa.info$parms,
      sa.syn.exp = sa.info$syn.exp,
      sp.parms   = sp.info$parms,
      sp.syn.exp = sp.info$syn.exp))
  }

