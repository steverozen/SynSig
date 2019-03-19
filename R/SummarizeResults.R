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


#' Assess/evaluate results from SigProfiler or SignatureAnalyzer
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
#' @param extracted.exp.path Path to extracted exposures file.
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
#' @export
#'
#' @importFrom utils capture.output sessionInfo

SummarizeSigOneSubdir <-
  function(third.level.dir,
           ground.truth.exposure.name,
           extracted.sigs.path,
           extracted.exp.path = NULL,
           # TODO(Steve): copy this to the summary and do analysis on how much
           # extracted signature contributs to exposures.

           read.extracted.sigs.fn,
           read.ground.truth.sigs.fn,
           write.cat.fn,
           plot.pdf.fn,
           overwrite = FALSE) {

    ## Output path - path to dump the ReadAndAnalyzeSigs() results
    outputPath <- paste0(third.level.dir, "/summary")

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

    if (FALSE) {
    copy.from <- paste0(third.level.dir,"/../ground.truth.syn.exposures.csv")
    if (!file.exists(copy.from)) {
      warning("Cannot find", copy.from, "\n\nSkipping\n\n")
    } else {
      copy.res <- file.copy(
        from = copy.from,
        to = paste0(third.level.dir,"/summary/"),
        overwrite = TRUE)
      if (!copy.res) cat("Copy from", copy.from, "to somewhere failed\n\n")
    }}

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

    capture.output(Sys.time(), sessionInfo(),
                   file = paste0(outputPath,"/log.txt"))

    invisible(sigAnalysis) # So we have something to check in tests
  }

#' Assess/evaluate results from SigProfiler-python (a.k.a. sigproextractor)
#'
#' @param third.level.dir Lowest level path to results, e.g.
#' \code{<top.dir>}\code{/sa.sa.96/sp.results/} or
#' \code{<top.dir>}\code{/sa.sa.96/sa.results/}
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
#' @param write.png If TRUE create png plots of the signatures.
#'
#'@param overwrite If TRUE overwrite existing directories and files.
#'
#' @export
#'
#' @importFrom ICAMS WriteCatSNS96 ReadCatSNS96
#' @importFrom utils capture.output sessionInfo
#' @importFrom grDevices png dev.off
#' @importFrom graphics par
#'
SummarizeSigOneSPSubdir <-
  function(third.level.dir,
           ground.truth.exposure.name = "ground.truth.syn.exposures.csv",
           write.png = FALSE,
           overwrite = FALSE) {

    # Location of SigProfiler output, which is our input
    # inputPath may change if sigproextractor updates!
    inputPath <- paste0(third.level.dir,"/SBS96/Suggested_Solution/De_Novo_Solution")
    stopifnot(dir.exists(inputPath))

    # Read in extracted signatures in sigproextractor txt format,
    # and convert it to ICAMS csv format.
    # Need special function to read in extracted signatures
    # Converted signatures will be included in the /summary folder.
    extractedSigs <- ReadSigProfilerSig96(paste0(inputPath,"/De_Novo_Solution_Signatures.txt"))
    extracted.sigs.path <- paste0(inputPath,"/extracted.signatures.PCAWG.format.csv")
    ICAMS::WriteCatSNS96(extractedSigs, extracted.sigs.path)

    retval <-
      SummarizeSigOneSubdir(
        third.level.dir = third.level.dir,
        ground.truth.exposure.name = ground.truth.exposure.name,
        extracted.sigs.path = extracted.sigs.path,
        read.extracted.sigs.fn = ReadCatSNS96,
        read.ground.truth.sigs.fn = ReadCatSNS96,
        write.cat.fn = WriteCatSNS96,
        plot.pdf.fn = PlotCatSNS96ToPdf,
        overwrite = overwrite)

    # Copy stability.pdf and result_stat.csv
    # generated by SigProfiler to summary/ folder
    # file.copy will return an "okay" flag, which equals to be TRUE if properly executed.
    # This is annoying, and I'll prevent this flag from printing it out
    copy.from.files <- paste0(third.level.dir,
                              c("/SBS96/All_Solution_Layer/L1/RE_vs_stabiliy_plot.pdf",
                                "/SBS96/All_Solution_Layer/L1/results_stat.csv"))
    for(copy.from in copy.from.files) {
      if (!file.exists(copy.from)) {
        warning("Cannot find", copy.from, "\n\nSkipping\n\n")
      } else {
        file.copy(# from = paste0(third.level.dir,"/SBS96/All_Solution_Layer/L1/stability.pdf"),
          from = copy.from,
          to = paste0(third.level.dir,"/summary/"),
          overwrite = TRUE)
      }
    }

    invisible(retval) # So we can test without looking at a file.
}

#' Summarize SigProfiler results in the sa.sa.96 and/or sp.sp subdirectories.
#'
#' @param top.dir The top directory of a conventional data structure containing
#' at least one of the subdirectories: sa.sa.96/sp.results and sp.sp/sp.results;
#' see further documentation elsewhere.
#'
#' @param sub.dir The subdirectory under \code{top.dir}, and containing a folder
#' named sp.results. By default, it contains both \code{c("sa.sa","sp.sp")}.
#' But you should specify \code{sub.dir = "sp.sp"} for \code{top.dir} with only
#' the \code{sp.sp} subdirectory
#' (as is the case for the correlated SBS1-and-SBS5-containing data sets).
#'
#' @param write.png If TRUE create png plots of the signatures.
#'
#' @export
#'
#' @details Results are put in standardized subdirectories of \code{top.dir}.

SummarizeSigProfiler <-
  function(top.dir, sub.dir = c("sa.sa.96","sp.sp"), write.png = FALSE) {

  ## If sub.dir are unexpected, throw an error
  expected.sub.dir <- c("sa.sa.96","sp.sp")
  if( !all(sub.dir %in% expected.sub.dir) ){ ## There are other sub-dirs than sa.sa.96 and sp.sp
    stop("sub.dir can only be one or two of c(\"sa.sa\",\"sp.sp\")!\n")
  }

  if("sa.sa.96" %in% sub.dir) {
    SummarizeSigOneSPSubdir(
      third.level.dir = paste0(top.dir, "/sa.sa.96/sp.results"),
      ground.truth.exposure.name = "ground.truth.syn.exposures.csv",
      write.png = write.png)
  }

  if("sp.sp" %in% sub.dir) {
    SummarizeSigOneSPSubdir(
      third.level.dir = paste0(top.dir, "/sp.sp/sp.results"),
      ground.truth.exposure.name = "ground.truth.syn.exposures.csv",
      write.png = write.png)
  }
}

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
#' @export
#'
#' @importFrom ICAMS WriteCatSNS96 ReadCatSNS96
#' @importFrom utils capture.output sessionInfo
#' @importFrom grDevices png dev.off
#' @importFrom graphics par
#
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
        read.extracted.sigs.fn = ReadCatCOMPOSITE,
        read.ground.truth.sigs.fn = ReadCatCOMPOSITE,
        write.cat.fn = WriteCatCOMPOSITE,
        plot.pdf.fn =  Plot96PartOfComposite, # NA, # Does not exist for COMPOSITE # maybe Plot96PartOfComposite
        overwrite = overwrite)

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
#' @param write.png If TRUE create png plots of the signatures.
#'
#' @param overwrite If TRUE overwrite existing directories and files.
#'
#' @export
#'
#' @importFrom ICAMS WriteCatSNS96 ReadCatSNS96
#' @importFrom utils capture.output sessionInfo
#' @importFrom grDevices png dev.off
#' @importFrom graphics par
#
SummarizeSigOneSA96Subdir <-
  function(third.level.dir,
           ground.truth.exposure.name = "ground.truth.syn.exposures.csv",
           which.run = "/best.run/",
           write.png = FALSE,
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
        read.extracted.sigs.fn = ReadCatSNS96,
        read.ground.truth.sigs.fn = ReadCatSNS96,
        write.cat.fn = WriteCatSNS96,
        plot.pdf.fn = PlotCatSNS96ToPdf,
        overwrite = overwrite)

    invisible(retval)
  }

#' Summarize all subdirectories of a major dataset.
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

#' Summarize all subdirectories of the correlated SBS1 / SBS5.
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

