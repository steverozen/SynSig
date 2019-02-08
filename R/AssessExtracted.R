#' Run \code{MatchSigs2Directions}, then
#' plot its results and write them as .csv files.
#'
#' @param ex.sigs Newly extracted signatures to be compared to gt.sigs
#            (actually, this is more general).
#
#' @param gt.sigs "Ground truth" signatures.
#'
#' @param exposure "Ground truth" exposures used generate the
#'   synthetic data from which \code{ex.sigs} were extracted.
#'
#' @return A list with the elements \code{avg}, \code{match1},
#' \code{match2} as for \code{SigSetSimilarity}, with \code{match1}
#'  being matches for the the extracted signatures (\code{ex.sigs})
#'  and \code{match2} being the
#'  matches for the ground truth signatures (\code{gt.sigs}). The return list
#'  also echos the input arguments \code{ex.sigs} and \code{gt.sigs}.
#
#' @export
#' @family signature matching functions

MatchSigsAndRelabel <-
  function(ex.sigs, gt.sigs, exposure) {

    if (is.null(colnames(ex.sigs))) {
      colnames(ex.sigs) <- paste0("Ex.", 1:ncol(ex.sigs))
    }

    # IMPORTANT Remove signatures that are not present in
    # the exposure from which the synthetic data were
    # generated
    exposed.sig.names <- rownames(exposure)[rowSums(exposure) > 0]
    # Make sure we do not have an signatures in exposures that
    # are not in gt.sigs.
    stopifnot(
      setequal(setdiff(exposed.sig.names, colnames(gt.sigs)), c()))
    gt.sigs <- gt.sigs[  , exposed.sig.names]

    sim <- MatchSigs2Directions(ex.sigs, gt.sigs)

    sim$extracted.with.no.best.match <-
      setdiff(colnames(ex.sigs), sim$match2$to)

    sim$ground.truth.with.no.best.match <-
      setdiff(colnames(gt.sigs), sim$match1$to)
    # TODO(Steve) Review documentation / explanation. Note that
    # e.g. SBS29 might have a best match (BI_COMPOSITE_SBS18_P)
    # but no BI signatures has SBS29 as its best match
    #

    # TODO(Steve) Document the complexity below; mostly it deals
    # with setting up plotting that is easy(?) to interpret.
    labels <- character(ncol(ex.sigs))
    names(labels) <- colnames(ex.sigs)
    nums <- NumFromId(sim$match1$to)
    reordered.ex <- colnames(ex.sigs)[order(nums)]
    ex.sigs.x <- ex.sigs[ , order(nums)]
    bestmatch.id <- sim$match1[reordered.ex, "to"]
    bestmatch.sim <- sim$match1[reordered.ex, "sim"]
    bestmatch.sim <- round(bestmatch.sim, digits=4)
    init.labels <-
      paste0(reordered.ex, " (", bestmatch.id, " ", bestmatch.sim, ")")
    names(init.labels) <- reordered.ex
    laggards <- setdiff(rownames(sim$match2), bestmatch.id)
    # Falling back to a loop here:
    for (lag in laggards) {
      my.ex.id  <- sim$match2[lag, "to"]
      my.ex.sim <- round(sim$match2[lag, "sim"], digits = 4)
      init.labels[my.ex.id] <-
        paste0(init.labels[my.ex.id],
               " (", lag, " ", my.ex.sim, ")")
    }
    colnames(ex.sigs.x) <- init.labels

    sim$ex.sigs <- ex.sigs.x
    sim$gt.sigs <- gt.sigs
    invisible(sim)
  }

#' @title Write exposure matrix to a file
#'
#' @param exposure.matrix Matrix of exposures
#'
#' @param file File to which to write the exposure matrix (as a CSV file)
#'
#' @export
#'
#' @importFrom utils write.csv
#'
WriteExposure <- function(exposure.matrix, file) {
  old.digits <- getOption("digits")
  options(digits = 22)
  write.csv(exposure.matrix, file, row.names = TRUE)
  options(digits = old.digits)
}

#' @title Read an exposure matrix to a file
#'
#' @param file CSV file containing an exposure matrix
#'
#' @return Matrix of exposures
#'
#' @export
#'
#' @importFrom utils read.csv
ReadExposure <- function(file) {
  return(read.csv(file, row.names = 1))
}

#' @title Assess how well extracted signatures match input signatures
#'
#' We assume that in many cases extraction programs will be run
#' outside of R on file inputs and will generate fill outputs.
#'
#' @param extracted.sigs File containing the extracted signature profiles.
#'
#' @param ground.truth.sigs File containing signature profiles from which the
#'  synthetic data were generated.
#'
#' @param read.extracted.sigs.fn Function to read the extracted signatures
#' into the appropriate standard internal representation. If NULL then
#' \code{read.ground.truth.sigs.fn} is used.
#'
#' @param read.ground.truth.sigs.fn Function to read the ground truth
#' signatures into the appropriate standard internal representation.
#'
#' @param ground.truth.exposures File containing the exposures from which
#'  the synthetic catalogs were generated.  This file is used to restrict
#'  assessment to only those signatures in \code{ground.truth.sigs}
#'  that were actually represented in the exposures.
#'
#' @return See \code{\link{MatchSigsAndRelabel}}
#'
#' @details Generates output files by calling
#' \code{\link{MatchSigsAndRelabel}}
#'
#' @export

ReadAndAnalyzeSigs <-
  function(extracted.sigs,
           ground.truth.sigs,
           ground.truth.exposures,
           read.extracted.sigs.fn = NULL,
           read.ground.truth.sigs.fn) {
  if (is.null(read.extracted.sigs.fn))
    read.extracted.sigs.fn <- read.ground.truth.sigs.fn

  ex.sigs <- read.extracted.sigs.fn(extracted.sigs)
  gt.sigs <- read.ground.truth.sigs.fn(ground.truth.sigs)
  exposure <- ReadExposure(ground.truth.exposures)
  # Rows are signatures, columns are samples.

  return(
    MatchSigsAndRelabel(ex.sigs, gt.sigs, exposure))
}

#' Create file names in a given directory
#'
#' The directory is provided by the global
#' variable \code{out.data.dir},
#' which \strong{must} be set by the user. If
#' \code{out.data.dir} is NULL then just return
#' \code{file.name}.
#'
#' @param file.name The name of the that will be
#' prefixed by \code{out.data.dir}.
#'
#' @return \code{file.name} prefixed by \code{out.data.dir}.
#'
#' @export
#'
OutDir <- function(file.name) {
  if (is.null(out.data.dir)) return(file.name)
  if (!dir.exists(out.data.dir)) {
      dir.create(out.data.dir)
  }
  return(paste0(out.data.dir, file.name))
}

#' Generate syntheic exposures from real exposures.
#'
#' Checkpoints the parameters and the synthetic
#' exposures to files. It also checks that the parameters
#' inferred from the synthetic data approximate those
#' inferred from \code{real.exp}.
#'
#' @param real.exp The actually exposures upon which to base
#' the parameters and synthetic exposures.
#'
#' @param num.syn.tumors Generate this number of synthetic tumors.
#'
#' @param file.prefix Prepend this to output filenames
#'  to indicate the organization of the data.
#'
#' @param sample.id.prefix Prefix for sample identifiers for the
#' synthetic samples.
#'
#' @return A list with elements:
#' \enumerate{
#'  \item \code{parms} The parameters inferred from \code{real.exp}.
#'  \item \code{syn.exp} The synthetic exposures generated from \code{parms}.
#' }
#'
#' @export
GenerateSynFromReal <-
  function(real.exp, num.syn.tumors, file.prefix, sample.id.prefix) {
    parms <- GetSynSigParamsFromExposures(real.exp)

    froot <- OutDir(file.prefix)

    WriteExposure(real.exp, paste0(froot, "real-exp.csv"))

    parm.file <- paste0(froot, "parms.csv")
    WriteSynSigParams(parms, parm.file)

    syn.exp <-
      GenerateSyntheticExposures(parms, num.syn.tumors, sample.id.prefix)

    WriteExposure(syn.exp, paste0(froot, "syn-exp.csv"))

    # Sanity check
    check.params <- GetSynSigParamsFromExposures(syn.exp)

    # sa.check.param should be similar to parms
    WriteSynSigParams(check.params, parm.file, append = TRUE)
    WriteSynSigParams(parms - check.params, parm.file,
                      append = TRUE)

    return(list(parms=parms, syn.exp=syn.exp))
  }
