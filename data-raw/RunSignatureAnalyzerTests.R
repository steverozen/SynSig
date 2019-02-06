## Run SignatureAnalyzer on sa.sa96 synthetic data.

# We need the SignatureAnalzer signatures, the exposures (so we know what
# signatures were used) and the synthetic data.

options( warn = 2 )

suppressMessages(library(data.table))
suppressMessages(library(lsa))
suppressMessages(library(SigMisc))
source("SignatureAnalyzer_interaction.R")
source("ICAMS_plotting.R")
source("ICAMS_set_cat_globals.R")
source("ICAMS_read_write_cat.R")
SetCatGlobals()
options( warn = 0 )
setwd("SignatureAnalzyer.052418/")
INPUT <- "INPUT_SignatureAnalyzer/"
source("SignatureAnalyzer.PCAWG.function.R")
setwd("..")
options( warn = 2 )

results.dir <- "SignatureAnalyzer_results_on_2018_12_31_v2_run5/"
in.data.dir <- "syn_test_data_2018_12_31_v2/"

ResultsFile <- function(file.name, ...) {
  if (is.null(results.dir)) return(paste0(file.name, ...))
  else {
    if (!dir.exists(results.dir)) {
      cat("Creating directory", results.dir, "to write file", file.name, "\n")
      dir.create(results.dir)
    }
    return(paste0(results.dir, file.name, ...))
  }
}

InFile <- function(file.name) {
  if (is.null(in.data.dir)) return(file.name)
  else {
    stopifnot(dir.exists(in.data.dir))
    return(paste0(in.data.dir, file.name))
  }
}

RunSignatureAnalyzerOnSyn <-
  function(syn.Rdata,  
           maxK = 30, tol = 1e-7, test.only = FALSE) {
    
    # Args
    #   syn.Rdata:   File to load, which contains the synthetic data
    #                in a single variable
    #   maxK:        Maximum number of signatures
    #   tol:         Tolerance for stopping SignatureAnalyzer
    #   test.only    If TRUE, then just use the first 10 columns of the
    #                synthetic data in syn.Rdata
    
    # Jaegil lowered tol to 1.e-05 for the PCAWG7 analysis.
   
    var <- load(InFile(syn.Rdata))
    stopifnot(length(var) == 1)
    syn.data <- eval(parse(text = var))
    
    if (test.only) syn.data <- syn.data[ , 1:10]
    
    info.root <- sub(".syn.Rdata", "",  syn.Rdata)
    dir <- paste0("tmp.", info.root)
    
    if (!dir.exists(dir)) {
      dir.create(dir) 
    }
    
    # TEMPORARY is a global required by SignatureAnalyzer
    TEMPORARY <<- paste0(dir, "/")
    
    # BayesNMF.L1W.L2H is defined by the statement
    # source("SignatureAnalyzer.PCAWG.function.R") above.
    out.data <- 
      BayesNMF.L1W.L2H(syn.data, 200000, 10, 5, tol, maxK, maxK, 1)
    
    sigs <- out.data[[1]]
    sigs <- sigs[   , colSums(sigs) > 0]
    
    save(sigs, file = ResultsFile(paste0(info.root, ".extracted.sigs.Rdata")))
    
    invisible(sigs)
  }


RunSignatureAnalyzerOnSyn("sa.sa.96.syn.Rdata")

RunSignatureAnalyzerOnSyn("sa.sa.COMPOSITE.syn.Rdata")

RunSignatureAnalyzerOnSyn("sa.sp.syn.Rdata")

RunSignatureAnalyzerOnSyn("sp.sp.syn.Rdata")

RunSignatureAnalyzerOnSyn("sp.sa.COMPOSITE.syn.Rdata")

# End of running SignatureAnalyzer; results should all be stored in .Rdata files.

AverageSimilarityOneTest <-
  function(extracted.sigs.to.load, gt.sigs, plot.fn, summary.file) {
    # Args:
    #   extracted.sigs.to.load: path to an .Rdata file containing
    #                           one variable contains the extracted
    #                           signatures from one test
    #   gt.sig:                 "ground truth" signatures from which
    #                           the test data were generated
    # Returns a list with the elments:
    #   avg, match1, match2: as for SigSetSimilarity, with sigs1 being the
    #                        extracted signatures and sigs2 being the 
    #                        ground truth signatures (gt.sig)
    #   
    #   ex.sigs: The extracted signatures stored in extracted.sigs.to.load
    #   
    #   gt.sig:  A copy of the input argument gt.sigs
    #   
    # Side effects:
    # 
    #   Plots the extracted signatures with identifiers showing the nearest 
    #   ground-truth signatures in the "results" directory, as defined by
    #   the global results.dir.

    info.root <- sub(".extracted.sigs.Rdata", "",  extracted.sigs.to.load)
    var <- load(ResultsFile(extracted.sigs.to.load))
    stopifnot(length(var) == 1)
    ex.sigs <- eval(parse(text = var))
    if (is.null(colnames(ex.sigs))) { 
      colnames(ex.sigs) <- paste0("Ex.", 1:ncol(ex.sigs))
    }
    
    sim <- MatchSigsThenWriteAndPlot(ex.sigs, gt.sigs, plot.fn, 
                                   file.prefix = 
                                     paste0(ResultsFile(info.root), "."))
    
    cat(info.root, sim$avg, ncol(gt.sigs), ncol(ex.sigs), 
        paste0("file://", in.data.dir, info.root, ".syn.Rdata\n"),
        sep = ",", file = summary.file, append = TRUE)
    return(sim)
  }

ReadCOMPOSITE <- function(file) {
  retval <- as.data.frame(fread(file))
  stopifnot(colnames(retval)[1] == "mutation.type")
  rownames(retval) <- retval[ , "mutation.type"]
  return(retval[ , -1])
}


summary.file <- ResultsFile("summary.csv")
cat("Test,AvgSimilarity,NumGTSigs,NumExtractedSigs,SynData\n",
    file=summary.file)
sa.sa.COMPOSITE.sigs.used <- 
  ReadCOMPOSITE(InFile("sa.sa.COMPOSITE.sigs.used.csv"))
sa.sa.COMPOSITE.results <-
  AverageSimilarityOneTest("sa.sa.COMPOSITE.extracted.sigs.Rdata",
                           sa.sa.COMPOSITE.sigs.used, 
                           Plot96PartOfComposite,
                           summary.file)

sa.sa.96.sigs.used <- ReadCat96(InFile("sa.sa.96.sigs.used.csv"), 
                                strict = TRUE)
sa.sa.96.results <-
  AverageSimilarityOneTest("sa.sa.96.extracted.sigs.Rdata", 
                           sa.sa.96.sigs.used,
                           Cat96ToPdf,
                           summary.file)

sa.sp.sigs.used <- ReadCat96(InFile("SP.sigs.used.in.sa.sp.csv"), 
                                strict = TRUE)
sa.sp.results <-
  AverageSimilarityOneTest("sa.sp.extracted.sigs.Rdata", 
                           sa.sp.sigs.used, 
                           Cat96ToPdf,
                           summary.file)

sp.sa.COMPOSITE.sigs.used <- 
  ReadCOMPOSITE(InFile("sp.sa.COMPOSITE.sigs.used.csv"))
sp.sa.COMPOSITE.results <-
  AverageSimilarityOneTest("sp.sa.COMPOSITE.extracted.sigs.Rdata",
                           sp.sa.COMPOSITE.sigs.used, 
                           Plot96PartOfComposite,
                           summary.file)

sp.sp.sigs.used <- ReadCat96(InFile("sp.sp.sigs.used.csv"), 
                             strict = TRUE)
sp.sp.results <- 
  AverageSimilarityOneTest("sp.sp.extracted.sigs.Rdata", 
                           sp.sp.sigs.used, 
                           Cat96ToPdf,
                           summary.file)
