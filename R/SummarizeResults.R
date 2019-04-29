CopyWithChecks <- function(from, to.dir, overwrite = FALSE) {
  if (!file.exists(from)) {
    warning("Cannot find", from, "\n\nSkipping\n\n")
  } else {
    copy.res <- file.copy(
      from = from, to = to.dir, overwrite = overwrite)
    if (!copy.res)
      cat("Copy from", from, "to directory", to.dir, "failed\n\n")
  }
}


#' Assess/evaluate results from Multiple software packages
#'
#' Note: Users should use sigproextractor(SigProfiler-Python) v0.0.5.43
#' and SignatureAnalyzer 2018-Apr-18
#'
#' @param third.level.dir Lowest level path to results, that is,
#' \code{top.dir}/sp.sp/sa.results/. or
#' \code{top.dir}/sa.sa.96/sa.results/
#' Here, \code{top.dir} refers to a top-level directory which contains the
#' full information of a synthetic dataset. (e.g. \code{syn.2.7a.7b.abst.v8})
#' This code depends on a conventional directory structure documented
#' elsewhere. However there should be a directory within the \code{third.level.dir}
#' which stores the software output.
#'
#' @param ground.truth.exposure.name File name which stores ground-truth exposures.
#' It defaults to be "ground.truth.syn.exposures.csv". This file can be found
#' in the \code{sub.dir}, i.e. \code{third.level.dir}/../
#'
#' @param extracted.sigs.path Path to extracted sigs file, e.g.
#' \code{<third.level.dir>/SBS96/Selected_Solution/De_Novo_Solution/signatures.PCAWG.format.csv}.
#'
#' @param attributed.exp.path Path to attributed exposures file.
#'
#' @param read.extracted.sigs.fn Function to read the extracted sigs file.
#' e.g. \code{ReadCatSNS96}
#'
#' @param read.ground.truth.sigs.fn Function to read the ground-truth sigs file.
#' e.g. \code{ReadCatSNS96}
#'
#' @param plot.pdf.fn If a function, use it to plot PDFs of the ground truth and
#' extracted signatures.
#'
#' @param write.cat.fn Function to write a catalog to disk, for example
#' [ICAMS]\code{WriteCatSNS96}
#'
#' @param overwrite If TRUE overwrite existing directories and files.
#'
#' @keywords internal
#'
#' @importFrom utils capture.output sessionInfo

SummarizeSigOneSubdir <-
  function(third.level.dir,
           ground.truth.exposure.name,
           extracted.sigs.path,
           attributed.exp.path = NULL,
           # TODO(Steve): copy this to the summary and do analysis on how much
           # extracted signature contributes to exposures.
           read.extracted.sigs.fn,
           read.ground.truth.sigs.fn,
           write.cat.fn,
           plot.pdf.fn,
           overwrite = FALSE) {

    ## Output path - path to dump the ReadAndAnalyzeSigs() results
    outputPath <- paste0(third.level.dir, "/summary")

    ## Analyze signature extraction similarity
    sigAnalysis <-
      ReadAndAnalyzeSigs(
        extracted.sigs = extracted.sigs.path,
        ground.truth.sigs =
          paste0(third.level.dir,"/../ground.truth.syn.sigs.csv"),
        ground.truth.exposures =
          paste0(third.level.dir,"/../", ground.truth.exposure.name),
        read.extracted.sigs.fn = read.ground.truth.sigs.fn,
        read.ground.truth.sigs.fn = read.ground.truth.sigs.fn)

    if (dir.exists(outputPath)) {
      if (!overwrite) stop(outputPath, " already exists")
    }
    suppressWarnings(dir.create(outputPath))

    # Copies ground.truth exposures from second.level.dir
    # to third.level.dir/summary.
    CopyWithChecks(
      from = paste0(third.level.dir,"/../ground.truth.syn.exposures.csv"),
      to.dir = paste0(third.level.dir,"/summary/"),
      overwrite = TRUE)

    # Writes bi-directional matching and cos.sim calculation
    write.csv(sigAnalysis$match1,
              file = paste(outputPath,"match1.csv",sep = "/"))
    write.csv(sigAnalysis$match2,
              file = paste(outputPath,"match2.csv",sep = "/"))

    # Writes ground truth and extracted signatures
    write.cat.fn(
      sigAnalysis$gt.sigs,
      path = paste(outputPath,"ground.truth.sigs.csv",sep = "/"))
    write.cat.fn(
      sigAnalysis$ex.sigs,
      path = paste(outputPath,"extracted.sigs.csv",sep = "/"))

    # Dumps other outputs into "other.results.txt"
    capture.output(
      cat("Average cosine similarity\n"),
      sigAnalysis$avg,
      cat("\nNumber of ground-truth signatures\n"),
      ncol(sigAnalysis$gt.sigs),
      cat("\nNumber of extracted signatures\n"),
      ncol(sigAnalysis$ex.sigs),
      cat("\nsigAnalysis$extracted.with.no.best.match\n"),
      sigAnalysis$extracted.with.no.best.match,
      cat("\nsigAnalysis$ground.truth.with.no.best.match\n"),
      sigAnalysis$ground.truth.with.no.best.match,
      file = paste0(outputPath,"/other.results.txt"))

    if (class(plot.pdf.fn) == "function") {
      # Output ground-truth sigs to a PDF file
      plot.pdf.fn(sigAnalysis$gt.sigs,
                  paste0(outputPath,"/ground.truth.sigs.pdf"),
                  type = "signature")

      # Output extracted sigs to a PDF file
      plot.pdf.fn(sigAnalysis$ex.sigs,
                  paste0(outputPath,"/extracted.sigs.pdf"),
                  type = "signature")
    }

    ## Analyze exposure attribution
    # To be compatible with PCAWG project which only studies
    # signature extraction not exposure attribution,
    # errors will not be thrown if exists("attributed.exp.path") == F.
    if(exists("attributed.exp.path")) {

      if(file.exists(attributed.exp.path)) {
        expDifference <- ReadAndAnalyzeExposures(
          extracted.sigs = extracted.sigs.path,
          ground.truth.sigs =
            paste0(third.level.dir,"/../ground.truth.syn.sigs.csv"),
          attributed.exp.path = attributed.exp.path,
          ground.truth.exposures =
            paste0(third.level.dir,"/../", ground.truth.exposure.name),
          read.extracted.sigs.fn = read.ground.truth.sigs.fn,
          read.ground.truth.sigs.fn = read.ground.truth.sigs.fn)

        # Write results of exposure attribution analysis
        write.csv(expDifference,
                  file = paste0(outputPath,"/exposureDifference.csv"),
                  quote = T)

        # Copy attributed exposures to summary folder.
        CopyWithChecks(attributed.exp.path,
                       paste0(outputPath,"/attributed.exposures.csv"),
                       overwrite = overwrite)
      }
      else {
        warning("Cannot find", attributed.exp.path, "\n\nSkipping\n\n")
      }
    }

    ## Log of system time and session info
    capture.output(Sys.time(), sessionInfo(),
                   file = paste0(outputPath,"/log.txt"))

    invisible(sigAnalysis) # So we have something to check in tests
  }





