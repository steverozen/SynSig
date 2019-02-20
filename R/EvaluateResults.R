#' Assess/evaluate results from SigProfiler or SignatureAnalyzer
#'
#' @param path Lowest level path to results, e.g.
#' ATopLevelDir/sa.sa.96/sp.results/.
#' This code depends on a conventional directory structure documented
#' elsewhere. However there should be a directory <path>/SBS96 which
#' stores SigProfiler results.
#'
#' @param software.name Name of the software from which the results are generated.
#' It can be the full name of the software (e.g. One in \code{c("sigproextractor","SigProfiler","SignatureAnalyzer")});
#' it can also be the acronym of the software (e.g. \code{c("sp","sa")} is short for
#' \code{c("SigProfiler","SignatureAnalyzer")}.
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
  function(path, software.name, ground.truth.exposures, ground.truth.sigs,
           write.png = FALSE) {

    ## Stop if the software name does not in the expected list
    software.correct.flag <- FALSE
    software.acronyms <- c("sp","sa")
    software.full.names <- list("sp" = c("sigproextractor","sigprofiler"),
                                "sa" = "signatureanalyzer")

    # If software.name has more than 1 parameter, stop!
    if(length(software.name) != 1) stop("You should provide only 1 software.name!\n")
    # If software.name is an acronym, then just convert it to lower case.
    if(tolower(software.name) %in% software.acronyms) {
      software.name <- tolower(software.name)
      software.correct.flag <- TRUE
    }
    # If software.name is a full name, convert it to its acronym.
    for(acro in names(software.full.names)) {
      if(tolower(software.name) %in% software.full.names[[acro]]){
        software.name <- acro
        software.correct.flag <- TRUE
      }
    }
    # If software.name is neither an acronym or a full name, raise an error.
    if(software.correct.flag != TRUE) stop("software.name is incorrect!\n")

    ## Input path - Location of tool output, which is our input
    # Location of SigProfiler output
    if(software.name == "sp") {
      inputPath <- paste0(path,"/SBS96/Selected_Solution/De_Novo_Solution")
    }
    # Location of SignatureAnalyzer output
    if(software.name == "sa"){
      inputPath <- paste0(path, "/sa.best.run/")
    }

    ## Output path - path to dump the ReadAndAnalyzeSigs() results
    outputPath <- paste0(path, "/evaluation")

    # For sigprofiler, we need special function to read in extracted signatures
    if(software.name == "sp") {
      extractedSigs <- ReadSigProfilerSig96(paste0(inputPath,"/signatures.txt"))
      WriteCat96(ct = extractedSigs,paste0(inputPath,"/signatures.PCAWG.format.csv"))
    }

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

    ## Old realization of PDF plot - keep it for debugging
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

    ## Output logs
    # Session Info
    capture.output(sessionInfo(),
      file = paste0(outputPath,"/sessionInfo.log"))
    # Time to finish the analysis
    capture.output(Sys.time,
      file = paste0(outputPath,"/finish.time.log"))
  }

#' Summarize SigProfiler results in the sa.sa.96 and sp.sp subdirectories
#'
#' @param top.dir The top directory of a conventional data structure containing
#' at least one of the subdirectories: sa.sa.96/sp.results and sp.sp/sp.results;
#' see further documentation elsewhere.
#'
#' @param sub.dir The subdirectory under \code{top.dir}, and containing a folder
#' named sp.results. By default, it equals to \code{c("sa.sa.96","sp.sp")}. But
#' you should specify \code{sub.dir = "sp.sp"} for clock-like SBS1-and-SBS5-containing
#' datasets.
#'
#' @param write.png If TRUE create png plots of the signatures.
#'
#' @export
#'
#' @details Results are put in standardized subdirectories of \code{top.dir}.

SummarizeSigProfiler <- function (top.dir, sub.dir = c("sa.sa.96","sp.sp"), write.png = FALSE) {

  if(0){
    stop("sub.dir can only be one or two of c(\"sa.sa\",\"sp.sp\")!\n")
  }

  if("sa.sa.96" %in% sub.dir) {
    SummarizeSigOneSubdir(
      paste0(top.dir, "/sa.sa.96/sp.results"),
      "sa.expsoure.csv",
      software.name = "sp",
      sa.96.sigs,
      write.png = write.png)
  }

  if("sp.sp" %in% sub.dir) {
    SummarizeSigOneSubdir(
      paste0(top.dir, "/sp.sp/sp.results"),
      "sp.expsoure.csv",
      software.name = "sp",
      sp.sigs,
      write.png = write.png)
  }

}
