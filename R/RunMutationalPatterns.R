#' Install MutationalPatterns from Bioconductor,
#' and its dependent package, NMF.
#'
#' @keywords internal
InstallMutationalPatterns <- function(){
  message("Installing MutationalPatterns from Bioconductor...")

  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("MutationalPatterns")

  if ("NMF" %in% rownames(installed.packages()) == FALSE)
    install.packages("NMF")
}


#' Run MutationalPatterns attribution on a spectra catalog file
#' and known signatures.
#'
#' @param input.catalog File containing input spectra catalog.
#' Columns are samples (tumors), rows are mutation types.
#'
#' @param gt.sigs.file File containing input mutational signatures.
#' Columns are signatures, rows are mutation types.
#'
#' @param read.catalog.function Function to read a catalog
#' (can be spectra or signature catalog): it takes a file path as
#' its only argument and returning a catalog as a numeric matrix.
#'
#' @param out.dir Directory that will be created for the output;
#' abort if it already exits.  Log files will be in
#' \code{paste0(out.dir, "/tmp")}.
#'
#' @param seed Specify the pseudo-random seed number
#' used to run MutationalPatterns. Setting seed can make the
#' attribution of MutationalPatterns repeatable.
#' Default: 1.
#'
#' @param test.only If TRUE, only analyze the first 10 columns
#' read in from \code{input.catalog}.
#' Default: FALSE
#'
#' @param overwrite If TRUE, overwrite existing output.
#' Default: FALSE
#'
#' @return The attributed exposure of \code{MutationalPatterns}, invisibly.
#'
#' @details Creates several
#'  files in \code{paste0(out.dir, "/sa.output.rdata")}. These are
#'  TODO(Steve): list the files
#'
#' @importFrom utils capture.output
#'
#' @export
#'
RunMutationalPatternsAttributeOnly <-
  function(input.catalog,
           gt.sigs.file,
           read.catalog.function,
           out.dir,
           seed = 1,
           test.only = FALSE,
           overwrite = FALSE) {

    ## Install MutationalPatterns, if not found in library
    if ("MutationalPatterns" %in% rownames(installed.packages()) == FALSE)
      InstallMutationalPatterns()


    ## Set seed
    set.seed(seed)
    seedInUse <- .Random.seed  ## Save the seed used so that we can restore the pseudorandom series
    RNGInUse <- RNGkind() ## Save the random number generator (RNG) used


    ## Read in spectra data from input.catalog file
    ## spectra: spectra data.frame in ICAMS format
    spectra <- read.catalog.function(input.catalog,
                                     strict = FALSE)
    if (test.only) spectra <- spectra[ , 1:10]

    ## Read in ground-truth signature file
    gtSignatures <- read.catalog.function(gt.sigs.file)

    ## Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }


    ## Derive exposure count attribution results.
    ## WARNING: MutationalPatterns can only do exposure attribution
    ## using SBS96 spectra catalog and signature catalog!
    exposureObject <-
      MutationalPatterns::fit_to_signatures(mut_matrix = spectra,
                                            signatures = gtSignatures)
    ## exposure attributions (in mutation counts)
    exposureCounts <- t(exposureObject$contribution)
    colnames(exposureCounts) <- paste("MP",1:ncol(exposureCounts),sep=".")
    ## Write exposure counts in ICAMS and SynSig format.
    WriteExposure(exposureCounts,
                  paste0(out.dir,"/attributed.exposures.csv"))

    ## Save seeds and session information
    ## for better reproducibility
    capture.output(sessionInfo(), file = paste0(out.dir,"/sessionInfo.txt")) ## Save session info
    write(x = seedInUse, file = paste0(out.dir,"/seedInUse.txt")) ## Save seed in use to a text file
    write(x = RNGInUse, file = paste0(out.dir,"/RNGInUse.txt")) ## Save seed in use to a text file

    ## Return attributed exposures
    invisible(exposureCounts)
  }



#' Run MutationalPatterns extraction and attribution on a spectra catalog file
#'
#' WARNING: MutationalPatterns can only do exposure attribution
#' using SBS96 spectra catalog and signature catalog!
#'
#' @param input.catalog File containing input spectra catalog.
#' Columns are samples (tumors), rows are mutation types.
#'
#' @param read.catalog.function Function to read a catalog
#' (can be spectra or signature catalog): it takes a file path as
#' its only argument and returning a catalog as a numeric matrix.
#'
#' @param write.catalog.function Function to write a catalog.
#'
#' @param out.dir Directory that will be created for the output;
#' abort if it already exits.  Log files will be in
#' \code{paste0(out.dir, "/tmp")}.
#'
#' @param CPU.cores Number of CPUs to use in running
#' MutationalPatterns. For a server, 30 cores would be a good
#' choice; while for a PC, you may only choose 2-4 cores.
#' By default (CPU.cores = NULL), the CPU.cores would be equal
#' to \code{(parallel::detectCores())/2}, total number of CPUs
#' divided by 2.
#'
#' @param seed Specify the pseudo-random seed number
#' used to run MutationalPatterns. Setting seed can make the
#' attribution of MutationalPatterns repeatable.
#' Default: 1.
#'
#' @param K,K.range \code{K} is the precise value for
#' the number of signatures active in spectra (K).
#' Specify \code{K} if you know precisely how many signatures
#' are active in the \code{input.catalog}, which is the
#' \code{ICAMS}-formatted spectra file.
#'
#' \code{K.range} is A numeric vector \code{(K.min,K.max)}
#' of length 2 which tell MutationalPatterns to search the best
#' signature number active in spectra, K, in this range of Ks.
#' Specify \code{K.range} if you don't know how many signatures
#' are active in the \code{input.catalog}.
#'
#' WARNING: You must specify only one of \code{K} or \code{K.range}!
#'
#' Default: NULL
#'
#' @param test.only If TRUE, only analyze the first 10 columns
#' read in from \code{input.catalog}.
#' Default: FALSE
#'
#' @param overwrite If TRUE, overwrite existing output.
#' Default: FALSE
#'
#' @return The attributed exposure of \code{MutationalPatterns}, invisibly.
#'
#' @details Creates several
#'  files in \code{out.dir}. These are:
#'  TODO(Steve): list the files
#'
#'  TODO(Wuyang)
#'
#' @importFrom utils capture.output
#'
#' @export

RunMutationalPatterns <-
  function(input.catalog,
           read.catalog.function,
           write.catalog.function,
           out.dir,
           CPU.cores = NULL,
           seed = 1,
           K = NULL,
           K.range = NULL,
           test.only = FALSE,
           overwrite = FALSE) {

    ## Check whether ONLY ONE of K or K.range is specified.
    bool1 <- is.numeric(K) & is.null(K.range)
    bool2 <- is.null(K) & is.numeric(K.range) & length(K.range) == 2
    stopifnot(bool1 | bool2)

    ## Install MutationalPatterns, if not found in library
    if ("MutationalPatterns" %in% rownames(installed.packages()) == FALSE)
      InstallMutationalPatterns()


    ## Set seed
    set.seed(seed)
    seedInUse <- .Random.seed  ## Save the seed used so that we can restore the pseudorandom series
    RNGInUse <- RNGkind() ## Save the random number generator (RNG) used


    ## Read in spectra data from input.catalog file
    ## spectra: spectra data.frame in ICAMS format
    spectra <- read.catalog.function(input.catalog,
                                     strict = FALSE)
    if (test.only) spectra <- spectra[ , 1:10]

    ## Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }

    ## CPU.cores specifies number of CPU cores to use.
    ## CPU.cores will be capped at 30.
    ## If CPU.cores is not specified, CPU.cores will
    ## be equal to the minimum of 30 or (total cores)/2
    if(is.null(CPU.cores)){
      CPU.cores = min(30,(parallel::detectCores())/2)
    } else {
      stopifnot(is.numeric(CPU.cores))
      if(CPU.cores > 30) CPU.cores = 30
    }


    ## Run NMF using ICAMS-formatted spectra catalog
    ## Determine the best number of signatures (K.best).
    ## If K is provided, use K as the K.best.
    ## If K.range is provided, determine K.best by doing raw extraction.
    if(bool1){
      gof_nmf <- NMF::nmf(spectra,
                          rank = K.range,     ## Rank specifies number of signatures you want to assess
                          method = "brunet",  ## "brunet" is the default NMF method in NMF package.
                          nrun = CPU.cores,
                          seed = seed)

      K.best <- K
      print(paste0("Assuming there are ",K.best," signatures active in input spectra."))
    }
    if(bool2){
      gof_nmf <- NMF::nmf(spectra,
                          rank = K.range,     ## Rank specifies number of signatures you want to assess
                          method = "brunet",  ## "brunet" is the default NMF method in NMF package.
                          nrun = CPU.cores,
                          seed = seed)


      ## Choose the best signature number (K.best) active in the spectra
      ## catalog (input.catalog).
      ##
      ## According to paper "A flexible R package for nonnegative matrix factorization"
      ## (Gaujoux & Seoighe, 2010), the most common approach to choose number of
      ## signature (K, a.k.a. rank in this paper) is to choose the smallest K for which
      ## cophenetic correlation coefficient starts decreasing.
      for(current.K in K.range)
      {
        current.summary <- NMF::summary(gof_nmf$fit[[as.character(current.K)]])
        current.cophenetic.coefficient <- current.summary["cophenetic"]

        next.summary <- NMF::summary(gof_nmf$fit[[as.character(current.K+1)]])
        next.cophenetic.coefficient <- next.summary["cophenetic"]

        if(current.cophenetic.coefficient > next.cophenetic.coefficient)
          break
      }
      K.best <- current.K ## Choose K.best as the smallest current.K whose cophenetic
                          ## is greater than cophenetic from (current.K+1).
      print(paste0("The best number of signatures is found.",
                   "It equals to: ",K.best))
    }


    ## Generates a list contain
    sigs_nmf <- MutationalPatterns::extract_signatures(spectra,K.best,CPU.cores)
    ## names(sigs_nmf)
    ## [1] "signatures"    "contribution"  "reconstructed"
    sigsRaw <- sigs_nmf$signatures ## un-normalized signature matrix
    extractedSignatures <- t(t(sigsRaw) / colSums(sigsRaw))   ## normalize each signature's sum to 1
    ## Add signature names for signature matrix extractedSignatures
    colnames(extractedSignatures) <-
      paste("MP",1:ncol(extractedSignatures),sep=".")
    ## Output extracted signatures in Duke-NUS format
    write.catalog.function(extractedSignatures,
                           paste0(out.dir,"/extracted.signatures.csv"))


    ## Derive exposure count attribution results.
    ## WARNING: MutationalPatterns can only do exposure attribution
    ## using SBS96 spectra catalog and signature catalog!
    exposureObject <- MutationalPatterns::fit_to_signatures(mut_matrix = spectra,
                                                            signatures = extractedSignatures)
    ## exposure attributions (in mutation counts)
    exposureCounts <- t(exposureObject$contribution)
    colnames(exposureCounts) <- paste("MP",1:ncol(exposureCounts),sep=".")
    ## Write exposure counts in ICAMS and SynSig format.
    WriteExposure(exposureCounts,
                  paste0(out.dir,"/attributed.exposures.csv"))

    ## Save seeds and session information
    ## for better reproducibility
    capture.output(sessionInfo(), file = paste0(out.dir,"/sessionInfo.txt")) ## Save session info
    write(x = seedInUse, file = paste0(out.dir,"/seedInUse.txt")) ## Save seed in use to a text file
    write(x = RNGInUse, file = paste0(out.dir,"/RNGInUse.txt")) ## Save seed in use to a text file

    ## Return a list of signatures and exposures
    invisible(list("signature" = extractedSignatures,
                   "exposure" = exposureCounts))
  }
