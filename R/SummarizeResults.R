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
#' @param extracted.sigs.path Full path of extracted sigs file, e.g.
#' \code{<third.level.dir>/SBS96/Selected_Solution/De_Novo_Solution/signatures.PCAWG.format.csv}.
#'
#' @param read.extracted.sigs.fn Function to read the extracted sigs file.
#' e.g. \code{ReadCat96}
#'
#' @param read.ground.truth.sigs.fn Function to read the ground-truth sigs file.
#' e.g. \code{ReadCat96}
#'
#' @param plot.pdf.fn If not NULL, use this function to plot PDFs of the ground truth and
#' extracted signatures.
#'
#' @param plot.png.fn If not NULL, use this function to plot PNGs of the ground truth and
#' extracted signatures.
#'
#' @export
#'
#' @importFrom ICAMS WriteCat96 ReadCat96 PlotCat96
#' @importFrom utils capture.output sessionInfo
#' @importFrom grDevices png dev.off
#' @importFrom graphics par
#'
SummarizeSigOneSubdir <-
  function(third.level.dir,
           ground.truth.exposure.name,
           extracted.sigs.path,
           read.extracted.sigs.fn,
           read.ground.truth.sigs.fn,
           write.cat.fn,
           plot.pdf.fn,
           plot.png.fn) {

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

    dir.create(outputPath)

    write.csv(sigAnalysis$match1, file = paste(outputPath,"match1.csv",sep = "/"))
    write.csv(sigAnalysis$match2, file = paste(outputPath,"match2.csv",sep = "/"))

    write.cat.fn(sigAnalysis$gt.sigs, path = paste(outputPath,"ground.truth.sigs.csv",sep = "/"))
    write.cat.fn(sigAnalysis$ex.sigs, path = paste(outputPath,"extracted.sigs.csv",sep = "/"))

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

    if (!is.null(plot.pdf.fn)) {
      # Output ground-truth sigs to a PDF file
      plot.pdf.fn(sigAnalysis$gt.sigs,
                  paste0(outputPath,"/ground.truth.sigs.pdf"),
                  type = "signature")

      # Output extracted sigs to a PDF file
      plot.pdf.fn(sigAnalysis$ex.sigs,
                  paste0(outputPath,"/extracted.sigs.pdf"),
                  type = "signature")
    }

    if(!is.null(plot.png.fn)) {
      # Output ground truth and extracted sigs to PNG files

      PngSigs <- function(sigs, sub.path) {
        out.sub.path <- paste0(outputPath,"/", sub.path, "/")
        dir.create(out.sub.path)

        PngOneSig <-  function(x) {
          png(paste0(out.sub.path, x, ".png"))
          par("mfrow"=c(1,1))
          plot.png.fn(sigs[ ,x, drop = FALSE], type = "signature",
                      grid = FALSE, xlabels = FALSE, cex = 0.6,
                      upper = FALSE)
          dev.off()
        }
        tmp <- lapply(colnames(sigs), PngOneSig)
      }

      PngSigs(sigAnalysis$gt.sigs, "ground.truth.sigs")
      PngSigs(sigAnalysis$ex.sigs, "extracted.sigs")

      if (FALSE) {
        PngOneGtSig <-  function(x) {
          png(paste0(outputPath,"/ground.truth.sigs/",x,".png"))
          par("mfrow"=c(1,1))
          plot.png.fn(sigAnalysis$gt.sigs[ ,x, drop = FALSE],
                      type = "signature",
                      grid = FALSE, xlabels = FALSE, cex = 0.6,
                      upper = FALSE)
          dev.off()
        }
        dir.create(paste0(outputPath,"/ground.truth.sigs/"))
        tmp <- lapply(colnames(sigAnalysis$gt.sigs), PngOneGtSig)

        PngOneExSig <-  function(x) {
          png(paste0(outputPath,"/extracted.sigs/",x,".png"))
          par("mfrow"=c(1,1))
          PlotCat96(sigAnalysis$ex.sigs[ ,x, drop = FALSE],
                    type = "signature",
                    grid = FALSE, xlabels = FALSE, cex = 0.6,
                    upper = FALSE)
          dev.off()
        }
        dir.create(paste0(outputPath,"/extracted.sigs/"))
        tmp <- lapply(colnames(sigAnalysis$ex.sigs), PngOneExSig)
      }
    }

    ## Logs
    # Session Info
    capture.output(sessionInfo(),
                   file = paste0(outputPath,"/sessionInfo.log"))
    # Date and time to finish the analysis
    capture.output(Sys.time(),
                   file = paste0(outputPath,"/finish.time.log"))

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
#' @export
#'
#' @importFrom ICAMS WriteCat96 ReadCat96 PlotCat96
#' @importFrom utils capture.output sessionInfo
#' @importFrom grDevices png dev.off
#' @importFrom graphics par
#'
SummarizeSigOneSPSubdir <-
  function(third.level.dir,
           ground.truth.exposure.name = "ground.truth.syn.exposures.csv",
           write.png = FALSE) {

    # Location of SigProfiler output, which is our input
    # inputPath may change if sigproextractor updates!
    inputPath <- paste0(third.level.dir,"/SBS96/Selected_Solution/De_Novo_Solution")
    # for SA inputPath <- paste0(path, "/sa.best.run/")

    # Need special function to read in extracted signatures
    extractedSigs <- ReadSigProfilerSig96(paste0(inputPath,"/signatures.txt"))
    extracted.sigs.path <- paste0(inputPath,"/signatures.PCAWG.format.csv")
    WriteCat96(ct = extractedSigs, extracted.sigs.path)

    SummarizeSigOneSubdir(
      third.level.dir = third.level.dir,
      ground.truth.exposure.name = ground.truth.exposure.name,
      extracted.sigs.path = extracted.sigs.path,
      read.extracted.sigs.fn = ReadCat96,
      read.ground.truth.sigs.fn = ReadCat96,
      write.cat.fn = WriteCat96,
      plot.png.fn = ifelse(write.png, PlotCat96, NULL),
      plot.pdf.fn = Cat96ToPdf)

    # Copy stability.pdf generated by SigProfiler to summary/ folder
    # file.copy will return an "okay" flag, which equals to be TRUE if properly executed.
    # This is annoying, and I'll prevent this flag from printing it out
    file.copy(from = paste0(third.level.dir,"/SBS96/All_Solution_Layer/L1/stability.pdf"),
              to = paste0(third.level.dir,"/summary/"),
              overwrite = TRUE)

    invisible(NULL)
}

#' Summarize SigProfiler results in the sa.sa.96 and/or sp.sp subdirectories
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

SummarizeSigProfiler <- function (top.dir, sub.dir = c("sa.sa.96","sp.sp"), write.png = FALSE) {

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

if (FALSE) {
SummarizeSigOneSPSubdir <-
  function(third.level.dir,
           ground.truth.exposure.name = "ground.truth.syn.exposures.csv",
           write.png = FALSE) {

    # Location of SignatureAnalyzer output, which is our input
    # inputPath may change if SignatureAnalyzer updates!
    inputPath <- paste0(third.level.dir,"/sa.best.run")

    # Need special function to read in extracted signatures
    extractedSigs <- ReadSigProfilerSig96(paste0(inputPath,"/signatures.txt"))
    extracted.sigs.path <- paste0(inputPath,"/signatures.PCAWG.format.csv")
    WriteCat96(ct = extractedSigs, extracted.sigs.path)

    SummarizeSigOneSubdir(
      third.level.dir = third.level.dir,
      ground.truth.exposure.name = ground.truth.exposure.name,
      extracted.sigs.path = extracted.sigs.path,
      read.extracted.sigs.fn = ReadCat96,
      read.ground.truth.sigs.fn = ReadCat96,
      write.cat.fn = WriteCat96,
      write.png = write.png)
  }
}

#' Summarize COMPOSITE results from SignatureAnalyzer
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
#' @export
#'
#' @importFrom ICAMS WriteCat96 ReadCat96 PlotCat96
#' @importFrom utils capture.output sessionInfo
#' @importFrom grDevices png dev.off
#' @importFrom graphics par
#
SummarizeSigOneSACOMPOSITESubdir <-
  function(third.level.dir,
           ground.truth.exposure.name = "ground.truth.syn.exposures.csv",
           which.run = "/sa.best.run/") {
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
        plot.pdf.fn = NULL, # Does not exist for COMPOSITE # maybe Plot96PartOfComposite
        plot.png.fn = NULL) # Does not exist for COMPOSITE

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
#' @export
#'
#' @importFrom ICAMS WriteCat96 ReadCat96 PlotCat96
#' @importFrom utils capture.output sessionInfo
#' @importFrom grDevices png dev.off
#' @importFrom graphics par
#
SummarizeSigOneSA96Subdir <-
  function(third.level.dir,
           ground.truth.exposure.name = "ground.truth.syn.exposures.csv",
           which.run = "/sa.best.run/",
           write.png = FALSE) {
    # Location of SigProfiler output, which is our input
    # inputPath may change if sigproextractor updates!
    inputPath <- paste0(third.level.dir, which.run)
    stopifnot(dir.exists(inputPath))

    retval <-
      SummarizeSigOneSubdir(
        third.level.dir = third.level.dir,
        ground.truth.exposure.name = ground.truth.exposure.name,
        extracted.sigs.path = paste0(inputPath,"/sa.output.sigs.csv"),
        read.extracted.sigs.fn = ReadCat96,
        read.ground.truth.sigs.fn = ReadCat96,
        write.cat.fn = WriteCat96,
        plot.pdf.fn = Cat96ToPdf,
        plot.png.fn = NULL) # Temporary

    invisible(retval)
  }

