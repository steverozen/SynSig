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
    ex.sigs.x <- ex.sigs[ , order(nums),drop = FALSE]
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

#' @title Assess how well extracted signatures match input signatures
#'
#' We assume that in many cases extraction programs will be run
#' outside of R on file inputs and will generate fill outputs.
#'
#' @param extracted.sigs Path to file containing the extracted signature profiles.
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


#' @title Assess how well attributed exposures match input exposures
#'
#' We assume that in many cases attribution programs will be run
#' outside of R on file inputs and will generate fill outputs.
#'
#' @param extracted.sigs Path to file containing the extracted signature profiles.
#'
#' @param ground.truth.sigs File containing signature profiles from which the
#'  synthetic data were generated.
#'
#' @param attributed.exp.path File containing mutation counts (exposures)
#' of synthetic tumors which are attributed to extracted or input signatures.
#'
#' @param ground.truth.exposures File containing the exposures from which
#'  the synthetic catalogs were generated.  This file is used to restrict
#'  assessment of signature exposures to only those signatures in
#'  \code{ground.truth.sigs} that were actually represented in the exposures.
#'
#' @param read.extracted.sigs.fn Function to read the extracted signatures
#' into the appropriate standard internal representation. If NULL then
#' \code{read.ground.truth.sigs.fn} is used.
#'
#' @param read.ground.truth.sigs.fn Function to read the ground truth
#' signatures into the appropriate standard internal representation.
#'
#' @return A \code{\link{data.frame}} recording:
#'
#' \code{Ground.truth.exposure}: sum of ground truth exposures of
#' all tumors to all ground-truth signatures.
#'
#' \code{Attributed.exposure}: sum of attributed exposures of
#' all tumors to all ground-truth signatures.
#' Here, attributed exposure of a tumor to a ground-truth
#' signature equals to the sum of the exposures of this tumor
#' to all extracted signatures which are most similar to
#' a ground-truth signature.
#' If there is no extracted signature resembling an ground-truth
#' signature, the attributed exposure of this ground-truth
#' signature will be \code{0}.
#'
#' \code{Absolute.difference}: sum of absolute difference between
#' ground-truth exposure and attributed exposure of all tumors
#' to all ground-truth signatures.
#'
#'
#' @details Generates output files by calling
#' \code{\link{MatchSigsAndRelabel}}
#'
#' @export

ReadAndAnalyzeExposures <-
  function(extracted.sigs,
           ground.truth.sigs,
           attributed.exp.path,
           ground.truth.exposures,
           read.extracted.sigs.fn = NULL,
           read.ground.truth.sigs.fn) {

  ## Bilaterally matching between ground-truth and extracted signatures
  sigMatch <- ReadAndAnalyzeSigs(extracted.sigs,
              ground.truth.sigs,
              ground.truth.exposures,
              read.extracted.sigs.fn = NULL,
              read.ground.truth.sigs.fn)


  ## Read in ground-truth and attributed exposures in ICAMS format
  gtExposures <- ReadExposure(ground.truth.exposures)
  attrExposures <- ReadExposure(attributed.exp.path)

  ## Names of ground-truth signatures
  gtSigsNames <- colnames(sigMatch$gt.sigs)

  ## Initialize an empty data.frame for exposure difference
  exposureSim <- data.frame(matrix(0,nrow = length(gtSigsNames),ncol = 3))
  rownames(exposureSim) <- gtSigsNames
  colnames(exposureSim) <- c("Ground.truth.exposure", ## Sum of all tumor's ground-truth exposure to gtSigsName
                             "Attributed.exposure", ## Sum of all tumor's attributed exposure to gtSigsName
                             "Absolute.difference") ## Sum of absolute difference of two exposure values for each tumor


  ## For each of the ground-truth signature, calculate the absolute difference
  ## between its input (ground-truth) exposure and its attributed exposure.
  ## Attributed exposure of a input signature equals to the sum of
  ## exposures of all extracted signatures which matches to
  ## this input signature.
  for(gtSigName in gtSigsNames){

    matchedExtrSigIndex <- which(sigMatch$match1[,1] == gtSigName)

    if(length(matchedExtrSigIndex) > 0) ## 1 or more extracted signatures match to gtSigName in match1
      matchedExtrSigName <- rownames(sigMatch$match1)[matchedExtrSigIndex]
    else ## No extracted signatures match to gtSigName
      matchedExtrSigName <- NULL

    for(index in 1:ncol(attrExposures)) { ## index refers to which tumor we are scrutinizing
      ## Each cycle traverses one tumor, and calculate the absolute difference
      ## between its attributed exposures and ground-truth exposures.
      gtExposureOneTumor <- gtExposures[gtSigName,index]
      attrExposureOneTumor <- ifelse(length(matchedExtrSigIndex) > 0,
                                     yes = sum(attrExposures[matchedExtrSigName,index]),
                                     no = 0)
      exposureSim[gtSigName,1] <- exposureSim[gtSigName,1] + gtExposureOneTumor
      exposureSim[gtSigName,2] <- exposureSim[gtSigName,2] + attrExposureOneTumor
      exposureSim[gtSigName,3] <- exposureSim[gtSigName,3] +
        abs(gtExposureOneTumor - attrExposureOneTumor)
    }

  }

  return(exposureSim)
}

