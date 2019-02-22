#' Assess/evaluate results from SigProfiler or SignatureAnalyzer
#'
#' @param third.level.dir Lowest level path to results, e.g.
#' ATopLevelDir/sa.sa.96/sp.results/.
#' This code depends on a conventional directory structure documented
#' elsewhere. However there should be a directory <third.level.dir>/SBS96 which
#' stores SigProfiler results.
#'
#' @param ground.truth.exposure.name One of "sa.exposure.csv" or "sp.exposure.csv".
#'
#' @param extracted.sigs.path XXXXX
#'
#' @param read.extracted.sigs.fn XXXXXX
#'
#' @param read.ground.truth.sigs.fn XXXXXXX
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
SummarizeSigOneSubdir <-
  function(third.level.dir,
           ground.truth.exposure.name,
           extracted.sigs.path,
           read.extracted.sigs.fn,
           read.ground.truth.sigs.fn,
           write.png = FALSE) {

    ## Output path - path to dump the ReadAndAnalyzeSigs() results
    outputPath <- paste0(third.level.dir, "/evaluation")

    sigAnalysis <-
      ReadAndAnalyzeSigs(
        extracted.sigs = extracted.sigs.path,
        ground.truth.sigs =
          paste0(third.level.dir,"/../ground.truth.sigs.csv"),
        ground.truth.exposures =
          paste0(third.level.dir,"/../../", ground.truth.exposure.name),
        read.ground.truth.sigs.fn = read.ground.truth.sigs.fn)

    dir.create(outputPath)

    write.csv(sigAnalysis$match1, file = paste(outputPath,"match1.csv",sep = "/"))
    write.csv(sigAnalysis$match2, file = paste(outputPath,"match2.csv",sep = "/"))

    WriteCat96(sigAnalysis$gt.sigs, path = paste(outputPath,"ground.truth.sigs.csv",sep = "/"))
    WriteCat96(sigAnalysis$ex.sigs, path = paste(outputPath,"extracted.sigs.csv",sep = "/"))

    capture.output(
      cat("\nAverage cosine similarity\n"),
      sigAnalysis$avg,
      cat("\nsigAnalysis$extracted.with.no.best.match\n"),
      sigAnalysis$extracted.with.no.best.match,
      cat("\nsigAnalysis$ground.truth.with.no.best.match\n"),
      sigAnalysis$ground.truth.with.no.best.match,
      file = paste0(outputPath,"other.results.txt"))

    ### Plot the input signatures
    dir.create(paste0(outputPath,"/ground.truth.sigs/"))

    PngOneSig <-  function(x) {
      png(paste0(outputPath,"/extracted.sigs/",x,".png"))
      par("mfrow"=c(1,1))
      PlotCat96(sigAnalysis$ex.sigs[ ,x, drop = FALSE],
                type = "signature",
                grid = FALSE, xlabels = FALSE, cex = 0.6,
                upper = FALSE)
      dev.off()
    }

    ## Optional: Output ground-truth sigs to PNG files
    if(write.png) {
      tmp <- lapply(colnames(sigAnalysis$gt.sigs), PngOneSig)
    }

    ## Output ground-truth sigs to a PDF file
    Cat96ToPdf(sigAnalysis$gt.sigs,
               paste0(outputPath,"/ground.truth.sigs/ground.truth.sigs.pdf"),
               type = "signature")

    ### Plot the extracted signatures
    dir.create(paste0(outputPath,"/extracted.sigs/"))

    ## Optional: Output extracted sigs to PNG files
    if(write.png) {
      tmp <- lapply(colnames(sigAnalysis$ex.sigs), PngOneSig)
    }

    ## Output extracted sigs to a PDF file

    Cat96ToPdf(sigAnalysis$ex.sigs,
               paste0(outputPath,
                      "/ground.truth.sigs/extracted.sigs.pdf"),
               type = "signature")

    ### Session Info
    capture.output(sessionInfo(),
                   file = paste0(outputPath,"/sessionInfo.log"))
  }


SummarizeSigOneSPSubdir <-
  function(third.level.dir,
           ground.truth.exposure.name,
           write.png = FALSE) {

    # Location of SigProfiler output, which is our input
    inputPath <- paste0(third.level.dir,"/SBS96/Selected_Solution/De_Novo_Solution")
    # for SA inputPath <- paste0(path, "/sa.best.run/")

    # Need special function to read in extracted signatures
    extractedSigs <- ReadSigProfilerSig96(paste0(inputPath,"/signatures.txt"))
    extracted.sigs.path <- paste0(inputPath,"/signatures.PCAWG.format.csv")
    WriteCat96(ct = extractedSigs, extracted.sigs.path)

    SummarizeSigOneSubdir(
      third.level.dir = third.level.dir,
      ground.truth.exposure.name = ground.truth.exposure.name,
      extracted.sigs = extracted.sigs.path,
      read.extracted.sigs.fn = ReadCat96,
      read.ground.truth.sigs.fn = ReadCat96,
      write.png = write.png)
}

#' Summarize SigProfiler results in the sa.sa.96 and sp.sp subdirectories
#'
#' @param top.dir The top directory of a conventional data structure containing
#' subdirectories sa.sa.96/sp.results and sp.sp/sp.results; see further
#' documentation elsewhere.
#'
#' @param write.png If TRUE create png plots of the signatures.
#'
#' @export
#'
#' @details Results are put in standardized subdirectories of \code{top.dir}.

SummarizeSigProfiler <- function (top.dir, write.png = FALSE) {

  SummarizeSigOneSPSubdir(
    third.level.dir = paste0(top.dir, "/sa.sa.96/sp.results"),
    ground.truth.exposure.name = "sa.expsoure.csv",
    write.png = write.png)

  SummarizeSigOneSPSubdir(
    paste0(top.dir, "/sp.sp/sp.results"),
    "sp.expsoure.csv",
    sp.sigs,
    write.png = write.png)
  }
