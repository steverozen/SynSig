#' Assess/evaluate SBS96 results from packages which can do
#' BOTH extraction and attribution,
#' excluding SigProfiler-Python and SignatureAnalyzer.
#'
#' Packages including but not limited to:
#' HDP, MutationalPatterns, sigfit,
#' SigneR, SomaticSignatures.
#'
#' @param third.level.dir Lowest level path to results, e.g.
#' \code{<top.dir>}\code{/sa.sa.96/SomaticSignatures.results/}
#' Here, \code{<top.dir>} refers to a top-level directory which contains the
#' full information of a synthetic dataset. (e.g. \code{syn.2.7a.7b.abst.v8})
#' This code depends on a conventional directory structure documented
#' elsewhere. For packages which can do both extraction and attribution,
#' we expect two files, \code{extracted.signatures.csv}
#' and \code{attributed.exposures.csv} are in the folder.
#'
#' @param ground.truth.exposure.name File name which stores ground-truth exposures;
#' defaults to \code{"ground.truth.syn.exposures.csv"}.
#' This file can be found in the \code{sub.dir}, i.e. \code{<third.level.dir>/../}
#'
#' @param overwrite If TRUE overwrite existing directories and files.
#'
#' @export
#'
#' @importFrom ICAMS WriteCatSNS96 ReadCatSNS96
#' @importFrom utils capture.output sessionInfo
#' @importFrom grDevices dev.off
#' @importFrom graphics par
#'
SummarizeSigOneExtrAttr96Subdir <-
  function(third.level.dir,
           ground.truth.exposure.name = "ground.truth.syn.exposures.csv",
           overwrite = FALSE) {

    # Location of SigProfiler output, which is our input
    # inputPath may change if sigproextractor updates!
    inputPath <- third.level.dir
    stopifnot(dir.exists(inputPath))

    # Specify the path of extracted signatures in ICAMS csv format.
    extracted.sigs.path <- paste0(inputPath,"/extracted.signatures.csv")

    # Specify the path of attributed exposures in SynSig csv format.
    attributed.exp.path <- paste0(inputPath,"/attributed.exposures.csv")

    # SummarizeSigOneSubdir will generate a "/summary" folder
    # under third.level.dir. Summarized results are dumped into
    # this folder.
    retval <-
      SummarizeSigOneSubdir(
        third.level.dir = third.level.dir,
        ground.truth.exposure.name = ground.truth.exposure.name,
        extracted.sigs.path = extracted.sigs.path,
        attributed.exp.path = attributed.exp.path,
        read.extracted.sigs.fn = ICAMS::ReadCatSNS96,
        read.ground.truth.sigs.fn = ICAMS::ReadCatSNS96,
        write.cat.fn = ICAMS::WriteCatSNS96,
        plot.pdf.fn = ICAMS::PlotCatSNS96ToPdf,
        overwrite = overwrite)

    invisible(retval) # So we can test without looking at a file.
  }





#' Assess/evaluate SBS96 results from packages which can
#' ONLY do exposure attribution.
#'
#' Packages including but not limited to:
#' deconstructSigs, YAPSA.
#'
#' Here, we excluded SignatureEstimation. Although it is also
#' a package with only attribution, but it has two attribution
#' algorithms. Therefore the naming of the results are slightly
#' different from the other two packages.
#'
#' @param third.level.dir Lowest level path to results, e.g.
#' \code{<top.dir>}\code{/sa.sa.96/SomaticSignatures.results/}
#' Here, \code{<top.dir>} refers to a top-level directory which contains the
#' full information of a synthetic dataset. (e.g. \code{syn.2.7a.7b.abst.v8})
#' This code depends on a conventional directory structure documented
#' elsewhere. For packages which can do both extraction and attribution,
#' we expect two files, \code{ground.truth.signatures.csv}
#' and \code{attributed.exposures.csv} are in the folder.
#'
#' @param ground.truth.exposure.name File name which stores ground-truth exposures;
#' defaults to \code{"ground.truth.syn.exposures.csv"}.
#' This file can be found in the \code{sub.dir}, i.e. \code{<third.level.dir>/../}
#'
#' @param overwrite If TRUE overwrite existing directories and files.
#'
#' @export
#'
#' @importFrom ICAMS WriteCatSNS96 ReadCatSNS96
#' @importFrom utils capture.output sessionInfo
#' @importFrom grDevices dev.off
#' @importFrom graphics par
#'
SummarizeSigOneAttr96Subdir <-
  function(third.level.dir,
           ground.truth.exposure.name = "ground.truth.syn.exposures.csv",
           overwrite = FALSE) {

    # Location of SigProfiler output, which is our input
    # inputPath may change if sigproextractor updates!
    inputPath <- third.level.dir
    stopifnot(dir.exists(inputPath))

    # Specify the path of extracted signatures in ICAMS csv format.
    ground.truth.sigs.path <- paste0(inputPath,"/../ground.truth.syn.sigs.csv")

    # SummarizeSigOneSubdir will generate a "/summary" folder
    # under third.level.dir. Summarized results are dumped into
    # this folder.
    retval <-
      SummarizeSigOneSubdir(
        third.level.dir = third.level.dir,
        ground.truth.exposure.name = ground.truth.exposure.name,
        extracted.sigs.path = ground.truth.sigs.path,
        attributed.exp.path = paste0(inputPath,"/attributed.exposures.csv"),
        read.extracted.sigs.fn = ReadCatSNS96,
        read.ground.truth.sigs.fn = ReadCatSNS96,
        write.cat.fn = WriteCatSNS96,
        plot.pdf.fn = PlotCatSNS96ToPdf,
        overwrite = overwrite)

    invisible(retval) # So we can test without looking at a file.
  }