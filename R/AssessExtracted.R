#' Run \code{MatchSigs2Directions}, then
#' plot its results and write them as .csv files.
#'
#' @param ex.sigs Newly extracted signatures to be compared to gt.sigs
#            (actually, this is more general)
#
#' @param gt.sigs "Ground truth" signatures
#'
#' @param plot.fn Function to plot signatures. For "typical" 96-channel
#            signatures use \code{Cat96ToPdf0}
#
#' @param file.prefix A prefix string for output files, can be used to put all
#                output files in one directory, e.g. by setting file.prefix
#                to "my_directory/something"
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
#'
#' @details Uses \code{plot.fn} to
#' plot the extracted signatures with identifiers showing the
#' nearest ground-truth signatures in a location defined by
#' \code{file.prefix}. Saves the match tables to a
#' .csv file with suffixes \code{match.extracted.to.gt.csv}
#' and \code{match.gt.to.extracted.csv}. \strong{IMPORTANT:}
#' This function does not check that all the \code{gt.sigs}
#' were "exposed" in generating the synthetic data from
#' which \code{ex.sigs} were generated. TODO(Steve): Possibly
#' add this check.
#'
MatchSigsThenWriteAndPlot <-
  function(ex.sigs, gt.sigs, plot.fn, file.prefix = '') {

    if (is.null(colnames(ex.sigs))) {
      colnames(ex.sigs) <- paste0("Ex.", 1:ncol(ex.sigs))
    }
    sim <- MatchSigs2Directions(ex.sigs, gt.sigs)

    # Write the "match" tables as .csv files
    fwriteDataFrame(sim$match1,
                    paste0(file.prefix, "match.extracted.to.gt.csv"),
                    rowname.name = "from")
    fwriteDataFrame(sim$match2,
                    paste0(file.prefix, "match.gt.to.extracted.csv"),
                    rowname.name = "from")

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
    laggards <- setdiff(sim$match2$from, bestmatch.id)
    # Falling back to a loop here:
    for (lag in laggards) {
      my.ex.id  <- sim$match2[lag, "to"]
      my.ex.sim <- round(sim$match2[lag, "sim"], digits = 4)
      init.labels[my.ex.id] <-
        paste0(init.labels[my.ex.id],
               " (", lag, " ", my.ex.sim, ")")
    }
    colnames(ex.sigs.x) <- init.labels
    plot.file <- paste0(file.prefix, "extracted.sigs.pdf")
    # cat("Plotting to", plot.file, "\n")
    plot.fn(ex.sigs.x, plot.file, type = "signature")

    sim$ex.sigs <- ex.sigs
    sim$gt.sigs <- gt.sigs
    invisible(sim)
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
#' @param read.ground.truth.sigs.fn Funtion to read the ground truth
#' signatures into the appropriate standard internal reprsenetation.
#'
#' @param ground.truth.exposures File containing the exposures from which
#'  the synthetic catalogs were generated.  This file is used to restrict
#'  assessment to only those signatures in \code{ground.truth.sigs}
#'  that were actually represented in the exposures.
#'
#' @param plot.fn Function to plot signatures. For "typical" 96-channel
#            signatures use \code{Cat96ToPdf0}
#
#' @param file.prefix A prefix string for output files, can be used to put all
#                output files in one directory, e.g. by setting file.prefix
#                to "my_directory/something"
#'
#' @return See \code{\link{MatchSigsThenWriteAndPlot}}
#'
#' @details Generates output files by calling
#' \code{\link{MatchSigsThenWriteAndPlot}}
#'
#' @export

ReadAndAnalyzeSigs <-
  function(extracted.sigs,
           ground.truth.sigs,
           ground.truth.exposures,
           read.extracted.sigs.fn = NULL,
           read.ground.truth.sigs.fn,
           plot.fn,
           file.prefix
           ) {
  if (is.null(read.extracted.sigs.fn))
    read.extracted.sigs.fn <- read.ground.truth.sigs.fn

  ex.sigs <- read.extracted.sigs.fn(extracted.sigs)
  gt.sigs <- read.ground.truth.sigs.fn(ground.truth.sigs)
  # TODO(Steve): IMPORTANT Read Exposures and remove signatures
  # that are not present in the exposures from the
  # ground.truth signatures.

  MatchSigsThenWriteAndPlot(ex.sigs, gt.sigs,
                            plot.fn, file.prefix)

}
