#' Summarize COMPOSITE results from SignatureAnalyzer.
#'
#' @param third.level.dir Lowest level path to results, that is
#' \code{<top.dir>}\code{/sa.sa.96/sa.results/},
#' \code{<top.dir>}\code{/sp.sp/sa.results/},
#'\code{<top.dir>}\code{/sa.sa.COMPOSITE/sa.results/}, or
#' \code{<top.dir>}\code{/sp.sa.COMPOSITE/sa.results/}.
#' Here, \code{<top.dir>} refers to a top-level directory which contains the
#' full information of a synthetic dataset. (e.g. \code{syn.2.7a.7b.abst.v8})
#' This code depends on a conventional directory structure documented
#' elsewhere. However there should be a directory
#' \code{<third.level.dir>}\code{/SBS96} which
#' stores SigProfiler results.
#'
#' @param ground.truth.exposure.name File name which stores ground-truth exposures;
#' defaults to \code{"ground.truth.syn.exposures.csv"}.
#' This file can be found in the \code{sub.dir}, i.e. \code{<third.level.dir>/../}
#'
#' @param which.run Name of subdirectory containing the run to summarize.
#'
#' @param overwrite If TRUE overwrite existing directories and files.
#'
#' @keywords internal
#'
#' @importFrom ICAMS WriteCatSNS96 ReadCatSNS96
#' @importFrom utils capture.output sessionInfo

SummarizeSigOneSACOMPOSITESubdir <-
  function(third.level.dir,
           ground.truth.exposure.name = "ground.truth.syn.exposures.csv",
           which.run = "/best.run/",
           overwrite = FALSE) {
    # Location of SigProfiler output, which is our input
    # inputPath may change if sigproextractor updates!
    inputPath <- paste0(third.level.dir, which.run)
    stopifnot(dir.exists(inputPath))

    retval <-
      SummarizeSigOneSubdir(
        third.level.dir = third.level.dir,
        ground.truth.exposure.name = ground.truth.exposure.name,
        extracted.sigs.path = paste0(inputPath,"/sa.output.sigs.csv"),
        attributed.exp.path = paste0(inputPath,"/sa.output.sigs.csv"),
        read.extracted.sigs.fn = ReadCatCOMPOSITE,
        read.ground.truth.sigs.fn = ReadCatCOMPOSITE,
        write.cat.fn = WriteCatCOMPOSITE,
        plot.pdf.fn =  Plot96PartOfComposite, # NA, # Does not exist for COMPOSITE # maybe Plot96PartOfComposite
        overwrite = overwrite)

    # TODO(Wuyang): Read in and Analyze attributed exposures,
    # if available.
    if(0){
      invisible(NULL)
      # attributed.exp.path
    }

    invisible(retval)
  }


#' Summarize 96-channel results from SignatureAnalyzer
#'
#' @param third.level.dir Lowest level path to results, that is
#' \code{<top.dir>}\code{/sa.sa.96/sa.results/},
#' \code{<top.dir>}\code{/sp.sp/sa.results/},
#'\code{<top.dir>}\code{/sa.sa.COMPOSITE/sa.results/}, or
#' \code{<top.dir>}\code{/sp.sa.COMPOSITE/sa.results/}.
#' Here, \code{<top.dir>} refers to a top-level directory which contains the
#' full information of a synthetic dataset. (e.g. \code{syn.2.7a.7b.abst.v8})
#' This code depends on a conventional directory structure documented
#' elsewhere. However there should be a directory
#' \code{<third.level.dir>}\code{/SBS96} which
#' stores SigProfiler results.
#'
#' @param ground.truth.exposure.name File name which stores ground-truth exposures;
#' defaults to \code{"ground.truth.syn.exposures.csv"}.
#' This file can be found in the \code{sub.dir}, i.e. \code{<third.level.dir>/../}
#'
#' @param which.run Name of subdirectory containing the run to summarize.
#'
#' @param overwrite If TRUE overwrite existing directories and files.
#'
#' @keywords internal
#'
#' @importFrom ICAMS WriteCatSNS96 ReadCatSNS96
#' @importFrom utils capture.output sessionInfo

SummarizeSigOneSA96Subdir <-
  function(third.level.dir,
           ground.truth.exposure.name = "ground.truth.syn.exposures.csv",
           which.run = "/best.run/",
           overwrite = FALSE) {
    # Location of SigProfiler output, which is our input
    # inputPath may change if sigproextractor updates!
    inputPath <- paste0(third.level.dir, which.run)

    if (!dir.exists(inputPath)) stop(inputPath, "does not exist")

    retval <-
      SummarizeSigOneSubdir(
        third.level.dir = third.level.dir,
        ground.truth.exposure.name = ground.truth.exposure.name,
        extracted.sigs.path = paste0(inputPath,"/sa.output.sigs.csv"),
        attributed.exp.path = paste0(inputPath,"/sa.output.exp.csv"),
        read.extracted.sigs.fn = ReadCatSNS96,
        read.ground.truth.sigs.fn = ReadCatSNS96,
        write.cat.fn = WriteCatSNS96,
        plot.pdf.fn = PlotCatSNS96ToPdf,
        overwrite = overwrite)

    # TODO(Wuyang): Read in and Analyze attributed exposures,
    # if available.
    if(0){
      invisible(NULL)
      # attributed.exp.path
    }

    invisible(retval)
  }

#' Summarize all subdirectories of SignatureAnalyzer results
#' on a major dataset.
#'
#' @param top.level.dir Path to top level directory.
#'
#' @param overwrite If TRUE overwrite existing summary files.
#'
#' @export

SignatureAnalyzerSummarizeTopLevel <-
  function(top.level.dir, overwrite = FALSE) {
    stopifnot(dir.exists(top.level.dir))

    assign("last.warning", NULL, envir = baseenv())

    options(warn = 2) # Warnings treated as errors

    sa.sa.96.dir <- paste0(top.level.dir, "/sa.sa.96/sa.results")
    stopifnot(dir.exists(sa.sa.96.dir))
    sp.sp.dir <- paste0(top.level.dir, "/sp.sp/sa.results")
    stopifnot(dir.exists(sp.sp.dir))
    sa.sa.COMPOSITE.dir <-
      paste0(top.level.dir, "/sa.sa.COMPOSITE/sa.results")
    stopifnot(dir.exists(sa.sa.COMPOSITE.dir))
    sp.sa.COMPOSITE.dir <-
      paste0(top.level.dir, "/sp.sa.COMPOSITE/sa.results")
    stopifnot(dir.exists(sp.sa.COMPOSITE.dir))

    CopyBestSignatureAnalyzerResult(sa.sa.96.dir, overwrite = overwrite)
    CopyBestSignatureAnalyzerResult(sp.sp.dir, overwrite = overwrite)
    CopyBestSignatureAnalyzerResult(sa.sa.COMPOSITE.dir, overwrite = overwrite)
    CopyBestSignatureAnalyzerResult(sp.sa.COMPOSITE.dir, overwrite = overwrite)

    retval <-
      list(sa.sa.96 =
             SummarizeSigOneSA96Subdir(
               sa.sa.96.dir, overwrite = overwrite),
           sp.sp =
             SummarizeSigOneSA96Subdir(
               sp.sp.dir, overwrite = overwrite),
           sa.sa.COMPOSITE =
             SummarizeSigOneSACOMPOSITESubdir(
               sa.sa.COMPOSITE.dir, overwrite = overwrite),
           sp.sa.COMPOSITE =
             SummarizeSigOneSACOMPOSITESubdir(
               sp.sa.COMPOSITE.dir, overwrite = overwrite))

    capture.output(print(retval), file = paste0(top.level.dir, "/retval.txt"))
    invisible(retval)
  }

#' Summarize all subdirectories of Signatureanalyzer results
#' on the correlated SBS1 / SBS5.
#'
#' This is special-purpose function to summarize results
#' from one in-silico experiment that examines how well
#' signatures can be extracted from synthetic tumors with
#' correlated SBS1 and SBS5.
#'
#' @param top.level.dir Path to top level directory.
#'
#' @param overwrite If TRUE overwrite existing directories and files.
#'
#' @export
SignatureAnalyzerSummarizeSBS1SBS5 <-
  function(top.level.dir, overwrite = FALSE) {
    stopifnot(dir.exists(top.level.dir))

    assign("last.warning", NULL, envir = baseenv())

    options(warn = 2) # Warnings treated as errors

    subdirs <-
      c("S.0.1.Rsq.0.1", "S.0.1.Rsq.0.2", "S.0.1.Rsq.0.3", "S.0.1.Rsq.0.6",
        "S.0.5.Rsq.0.1", "S.0.5.Rsq.0.2", "S.0.5.Rsq.0.3", "S.0.5.Rsq.0.6",
        "S.1.Rsq.0.1",   "S.1.Rsq.0.2",   "S.1.Rsq.0.3",   "S.1.Rsq.0.6")

    Summarize1 <- function(pseudo.top.level) { # pseudo.third.level is e.g. S.0.1.Rsq.0.1
      sp.sp.path <- paste0(top.level.dir, "/", pseudo.top.level, "/sp.sp/")
      if (!dir.exists(sp.sp.path)) stop(sp.sp.path, "does not exist")
      sa.results.path <- paste0(sp.sp.path, "sa.results/")
      if (!dir.exists(sa.results.path)) stop(sa.results.path, "does not exist")
      CopyBestSignatureAnalyzerResult(sa.results.path, overwrite = overwrite)
      SummarizeSigOneSA96Subdir(sa.results.path, overwrite = overwrite)
    }

    retval <- lapply(subdirs, Summarize1)

    capture.output(print(retval), file = paste0(top.level.dir, "/retval.txt"))
    invisible(retval)
}
