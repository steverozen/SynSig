#' Assess/evaluate results from SigProfiler or SignatureAnalyzer
#'
#' @param path Lowest level path to results, e.g.
#' ATopLevelDir/sa.sa.96/sp.results/.
#' This code depends on a conventional directory structure documented
#' elsewhere. However there should be a directory <path>/SBS96 which
#' stores SigProfiler results.
#'
#' @param ground.truth.exposures One of "sa.exposure.csv" or "sp.exposure.csv".
#'
#' @param ground.truth.sigs One of \code{sp.sigs} or \code{sa.96.sigs}.
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
  function(path, ground.truth.exposures, ground.truth.sigs,
           write.png = FALSE) {

    # Location of SigProfiler output, which is our input
    inputPath <- paste0(path,"/SBS96/Selected_Solution/De_Novo_Solution")
    # for SA inputPath <- paste0(path, "/sa.best.run/")

    ## Output path - path to dump the ReadAndAnalyzeSigs() results
    outputPath <- paste0(path, "/evaluation")

    # Need special function to read in extracted signatures
    extractedSigs <- ReadSigProfilerSig96(paste0(inputPath,"/signatures.txt"))
    WriteCat96(ct = extractedSigs,paste0(inputPath,"/signatures.PCAWG.format.csv"))

    sigAnalysis <-
      ReadAndAnalyzeSigs(
        extracted.sigs =
          paste0(inputPath,"/signatures.PCAWG.format.csv"),
        ground.truth.sigs = ground.truth.sigs,
        ground.truth.exposures =
          paste0(path,"/../../", ground.truth.exposures),
        read.ground.truth.sigs.fn = ReadCat96)

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

    if (FALSE) {
    pdf(paste0(outputPath,"/ground.truth.sigs/ground.truth.sigs.pdf"))
    par(mfrow=c(8,1))
    tmp <-
      lapply(colnames(sigAnalysis$gt.sigs),
             function(x) {
               PlotCat96(sigAnalysis$gt.sigs[ ,x, drop = FALSE],
                         type = "signature",
                         grid = FALSE, xlabels = FALSE, cex = 0.6,
                         upper = FALSE)
             })
    dev.off()
    }

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


    if (FALSE) {
    pdf(paste0(outputPath,"/extracted.sigs/extracted.sigs.pdf"))
    par(mfrow=c(8,1))
    tmp <-
      lapply(colnames(sigAnalysis$ex.sigs),
             function(x) {
               PlotCat96(sigAnalysis$ex.sigs[ ,x, drop = FALSE],
                         type = "signature",
                         grid = FALSE, xlabels = FALSE, cex = 0.6,
                         upper = FALSE)
             })
    dev.off()
    }

    ### Session Info
    capture.output(sessionInfo(),
                   file = paste0(outputPath,"/sessionInfo.log"))
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

  SummarizeSigOneSubdir(
    paste0(top.dir, "/sa.sa.96/sp.results"),
    "sa.expsoure.csv",
    sa.96.sigs,
    write.png = write.png)

  SummarizeSigOneSubdir(
    paste0(top.dir, "/sp.sp/sp.results"),
    "sp.expsoure.csv",
    sp.sigs,
    write.png = write.png)
  }
