#' Install signeR from Bioconductor
#'
#' @keywords internal
InstallsigneR <- function(){
  message("Installing signeR from Bioconductor...\n")

  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("signeR")
}



#' Run signeR extraction and attribution on a spectra catalog file
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
#' @param seedNumber Specify the pseudo-random seed number
#' used to run signeR. Setting seed can make the
#' attribution of signeR repeatable.
#' Default: 1.
#'
#' @param K,K.range \code{K} is the precise value for
#' the number of signatures active in spectra (K).
#' Specify \code{K} if you know precisely how many signatures
#' are active in the \code{input.catalog}, which is the
#' \code{ICAMS}-formatted spectra file.
#'
#' \code{K.range} is A numeric vector \code{(K.min,K.max)}
#' of length 2 which tell signeR to search the best
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
#' @return The attributed exposure of \code{signeR}, invisibly.
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

RunsigneR <-
  function(input.catalog,
           read.catalog.function,
           write.catalog.function,
           out.dir,
           seedNumber = 1,
           K = NULL,
           K.range = NULL,
           test.only = FALSE,
           overwrite = FALSE) {

    ## Check whether ONLY ONE of K or K.range is specified.
    bool1 <- is.numeric(K) & is.null(K.range)
    bool2 <- is.null(K) & is.numeric(K.range) & length(K.range) == 2
    stopifnot(bool1 | bool2)

    ## Install signeR, if not found in library
    if ("signeR" %in% rownames(installed.packages()) == FALSE)
      InstallsigneR()


    ## Set seed
    set.seed(seedNumber)
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

    ## Convert ICAMS-formatted spectra and signatures into signeR format.
    ## The program requires the input to be an TRANSPOSED mutation count matrix.
    convSpectra <- spectra
    dimnames(convSpectra) <- dimnames(spectra)
    sample.number <- dim(spectra)[2]
    convSpectra <- t(convSpectra)

    ## Determine the best number of signatures (K.best).
    ## If K is provided, use K as the K.best.
    ## If K.range is provided, determine K.best by doing raw extraction.
    if(bool1){
      extractionObject <- signeR::signeR(M=convSpectra,  ## M: Mutation spectra you want to decompose
                                         #Opport = NULL, ## Opport: Abundance (Opportunity) matrix for the spectra (optional)
                                         nsig=K)         ## nsig: Number of signatures (K)
      K.best <- K
      print(paste0("Assuming there are ",K.best," signatures active in input spectra."))
    }
    if(bool2){
      ## Extraction and attribution when number of signatures (K) is not known:
      ## automatically determine best number of signatures,
      ## based on median Bayesian Information Criterion (BIC).
      ## Step 1: do raw extraction and attribution to test burn-in (1000 times)
      ## and test sampling (1000) for possible signature numbers (N);
      ## Step 2: Compare the BIC of these Ns and determine the best number of signatures (Nbest);
      ## Step 3: do precise extraction and attribution (burn-in: 10000, sampling: 2000);
      extractionObject <- signeR::signeR(M=convSpectra,  ## M: Mutation spectra you want to decompose
                                         #Opport = NULL, ## Opport: Abundance (Opportunity) matrix for the spectra (optional)
                                         nlim=K.range)   ## nlim: Minimal and maximal number of signatures (K.range)

      ## Record best number of signatures, and verify this choice using BIC-plot
      K.best <- extractionObject$Nsign
      print(paste0("The best number of signatures is found.",
                   "It equals to: ",K.best))
      pdf(paste0(out.dir,"/Nsig.BIC.plot.pdf"))
      signeR::BICboxplot(extractionObject)
      dev.off()
    }


    ## Output extracted signatures in Duke-NUS format
    extractedSignaturesRaw <- extractionObject$Phat
    colnames(extractedSignaturesRaw) <- paste("sig",seq(1,K.best),sep = ".")
    ## Normalize the extracted signatures so that frequencies of each signature sums up to 1
    extractedSignatures <- apply(extractedSignaturesRaw,2, function(x) x/sum(x))
    rownames(extractedSignatures) <- rownames(spectra)
    ## Write extracted signatures
    write.catalog.function(extractedSignatures,
                           paste0(out.dir,"/extracted.signatures.csv"))


    ## Derive exposure count attribution results.
    exposureCounts <- extractionObject$Ehat ## Unnormalized exposures
    rownames(exposureCounts) <- paste("sig",seq(1,K.best),sep = ".") ## Assign row names of exposure matrix as names of signatures
    colnames(exposureCounts) <- colnames(spectra) ## Assign column names of exposure matrix as names of tumors
    ## Normalize the attributed counts so that each column represents exposure of a signature
    for(ii in 1:ncol(exposureCounts)) {
      exposureCounts[,ii] <- exposureCounts[,ii] / sum(exposureCounts[,ii])
      exposureCounts[,ii] <- exposureCounts[,ii] * colSums(spectra)[ii]
    }
    ## Save exposure attribution results
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
