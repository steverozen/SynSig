#' Install sigfit from GitHub,
#' and its dependent package, rstan.
#'
#' @keywords internal
Installsigfit <- function(){
  message("Installing sigfit from GitHub kgori/sigfit ...\n")

  if ("devtools" %in% rownames(installed.packages()) == FALSE)
    install.packages("devtools")
  devtools::install_github("kgori/sigfit", args = "--preclean", build_vignettes = TRUE)

  ## Install package rstan, which is required in the run.
  if ("rstan" %in% rownames(installed.packages()) == FALSE)
    install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
}


#' Run sigfit attribution on a spectra catalog file
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
#' @param seedNumber Specify the pseudo-random seed number
#' used to run sigfit. Setting seed can make the
#' attribution of sigfit repeatable.
#' Default: 1.
#'
#' @param test.only If TRUE, only analyze the first 10 columns
#' read in from \code{input.catalog}.
#' Default: FALSE
#'
#' @param overwrite If TRUE, overwrite existing output.
#' Default: FALSE
#'
#' @return The attributed exposure of \code{sigfit}, invisibly.
#'
#' @details Creates several
#'  files in \code{paste0(out.dir, "/sa.output.rdata")}. These are
#'  TODO(Steve): list the files
#'
#' @importFrom utils capture.output
#'
#' @export
#'
RunsigfitAttributeOnly <-
  function(input.catalog,
           gt.sigs.file,
           read.catalog.function,
           out.dir,
           seedNumber = 1,
           test.only = FALSE,
           overwrite = FALSE) {

    ## Install MutationalPatterns, if not found in library
    if ("sigfit" %in% rownames(installed.packages()) == FALSE)
      Installsigfit()


    ## Set seed
    set.seed(seedNumber)
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

    ## Transpose spectra catalog and ground-truth
    ## signature catalog so that sigfit can analyse them.
    convSpectra <- t(spectra)
	convGtSigs <- t(gtSignatures)


    ## Derive exposure count attribution results.
    mcmc_samples_fit <- sigfit::fit_signatures(counts = convSpectra,
                                               signatures = convGtSigs,
                                               iter = 2000,
                                               warmup = 1000,
                                               chains = 1,
                                               seed = seedNumber)

    ## exposuresObj$mean contain the mean of attributed exposures across multiple
    ## MCMC samples. Note that attributed exposure in exposuresObj$mean are un-normalized.
    exposuresObj <- sigfit::retrieve_pars(mcmc_samples_fit,
                                          par = "exposures",
                                          hpd_prob = 0.90)


    ## mutation count of each tumor
    sum_mutation_count <- rowSums(convSpectra)
    ## Multiply relative exposure with total mutation count of each tumor
    exposureCounts <- exposuresObj$mean * sum_mutation_count ## i-th row will be multiplied by i-th element (row) of sum_mutation_count

    ## Change signature names for exposure matrix exposureCounts:
    ## E.g., replace "Signature A" with "Signature.A".
    colnames(exposureCounts) <-
      gsub(pattern = " ",replacement = ".",colnames(exposureCounts))



    ## Write exposure counts in ICAMS and SynSig format.
    WriteExposure(exposureCounts,
                  paste0(out.dir,"/attributed.exposures.csv"))

    ## Copy ground.truth.sigs to out.dir
    file.copy(from = gt.sigs.file,
              to = paste0(out.dir,"/ground.truth.signatures.csv"),
              overwrite = overwrite)

    ## Save seeds and session information
    ## for better reproducibility
    capture.output(sessionInfo(), file = paste0(out.dir,"/sessionInfo.txt")) ## Save session info
    write(x = seedInUse, file = paste0(out.dir,"/seedInUse.txt")) ## Save seed in use to a text file
    write(x = RNGInUse, file = paste0(out.dir,"/RNGInUse.txt")) ## Save seed in use to a text file

    ## Return attributed exposures
    invisible(exposureCounts)
  }



#' Run sigfit extraction and attribution on a spectra catalog file
#'
#' WARNING: sigfit can only do exposure attribution
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
#' sigfit. For a server, 30 cores would be a good
#' choice; while for a PC, you may only choose 2-4 cores.
#' By default (CPU.cores = NULL), the CPU.cores would be equal
#' to \code{(parallel::detectCores())/2}, total number of CPUs
#' divided by 2.
#'
#' @param seedNumber Specify the pseudo-random seed number
#' used to run sigfit. Setting seed can make the
#' attribution of sigfit repeatable.
#' Default: 1.
#'
#' @param K,K.range \code{K} is the precise value for
#' the number of signatures active in spectra (K).
#' Specify \code{K} if you know precisely how many signatures
#' are active in the \code{input.catalog}, which is the
#' \code{ICAMS}-formatted spectra file.
#'
#' \code{K.range} is A numeric vector \code{(K.min,K.max)}
#' of length 2 which tell sigfit to search the best
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
#' @return The attributed exposure of \code{sigfit}, invisibly.
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

Runsigfit <-
  function(input.catalog,
           read.catalog.function,
           write.catalog.function,
           out.dir,
           CPU.cores = NULL,
           seedNumber = 1,
           K = NULL,
           K.range = NULL,
           test.only = FALSE,
           overwrite = FALSE) {

    ## Check whether ONLY ONE of K or K.range is specified.
    bool1 <- is.numeric(K) & is.null(K.range)
    bool2 <- is.null(K) & is.numeric(K.range) & length(K.range) == 2
    stopifnot(bool1 | bool2)

    ## Install sigfit, if not found in library
    if ("sigfit" %in% rownames(installed.packages()) == FALSE)
      Installsigfit()


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

    ## For faster computation, enable parallel computing in rstan.
    ##
    ## Allow precompiled rstan program to be written temporarily,
    ## so that parallel call may be faster
    rstan::rstan_options(auto_write = TRUE)
    ## Use CPU.cores cores for parallel computing
    options(mc.cores = CPU.cores)

    ## Transpose spectra catalog so that sigfit
    ## can analyse it.
    convSpectra <- t(spectra)

    ## Determine the best number of signatures (K.best).
    ## If K is provided, use K as the K.best.
    ## If K.range is provided, determine K.best by doing raw extraction.
    if(bool1){

      K.best <- K
      print(paste0("Assuming there are ",K.best," signatures active in input spectra."))
    }
    if(bool2){

      ## Choose the best signature number (K.best) active in the spectra
      ## catalog (input.catalog).
      ## Raw extraction: estimate most likely number of signatures
      ## (Nsig.max + 1) number of elements in mcmc_samples_extr
      ## The first Nsig.min number of elements are NULL elements
      ## Nsig.min+1 to Nsig.max elements are list elements of two elements: $data and $result
      ## The last element is the best signature number
      mcmc_samples_extr <-
        sigfit::extract_signatures(counts = convSpectra,   ## The spectra matrix required in signature extraction
                                   nsignatures = K.range,  ## The possible number of signatures a spectra may have.
                                   model = "nmf",          ## Method to use: we choose "nmf" by default. We can also choose "emu"
                                   iter = 1000,            ## Number of iterations in the run
                                   seed = seedNumber)

      K.best <- mcmc_samples_extr$best ## Choose K.best
      print(paste0("The best number of signatures is found.",
                   "It equals to: ",K.best))

      ## Remove the raw extraction object,
      ## which is extremely large (~32G)
      rm(mcmc_samples_extr)
      gc() ## Do garbage collection to recycle RAM
    }


    ## Precise extraction:
    ## Specifying number of signatures, and iterating more times to get more precise extraction
    ## Return a list with two elements: $data and $result

    mcmc_samples_extr_precise <-
      sigfit::extract_signatures(counts = convSpectra,   ## The spectra matrix required in signature extraction
                                 nsignatures = K.best,   ## The possible number of signatures a spectra may have.
                                 model = "nmf",          ## Method to use: we choose "nmf" by default. We can also choose "emu"
                                 iter = 5000,            ## Number of iterations in the run
                                 seed = seedNumber)

    extrSigsObject <- sigfit::retrieve_pars(mcmc_samples_extr_precise,
                                            par = "signatures")
    ## mean extracted signatures
    extractedSignatures <- t(extrSigsObject$mean)

    ## Change signature names for signature matrix extractedSignatures:
    ## E.g., replace "Signature A" with "Signature.A".
    colnames(extractedSignatures) <-
      gsub(pattern = " ",replacement = ".",colnames(extractedSignatures))

    ## Write extracted signatures into a ICAMS signature catalog file.
    write.catalog.function(extractedSignatures,
                           paste0(out.dir,"/extracted.signatures.csv"))


    ## Derive exposure count attribution results.
    ## WARNING: sigfit can only do exposure attribution
    ## using SBS96 spectra catalog and signature catalog!
    mcmc_samples_fit <- sigfit::fit_signatures(counts = convSpectra,
                                               signatures = extrSigsObject$mean,
                                               iter = 2000,
                                               warmup = 1000,
                                               chains = 1,
                                               seed = 1)

    ## exposuresObj$mean contain the mean of attributed exposures across multiple
    ## MCMC samples. Note that attributed exposure in exposuresObj$mean are un-normalized.
    exposuresObj <- sigfit::retrieve_pars(mcmc_samples_fit,
                                          par = "exposures",
                                          hpd_prob = 0.90)


    ## mutation count of each tumor
    sum_mutation_count <- rowSums(convSpectra)
    ## Multiply relative exposure with total mutation count of each tumor
    exposureCounts <- exposuresObj$mean * sum_mutation_count ## i-th row will be multiplied by i-th element (row) of sum_mutation_count

    ## Change signature names for exposure matrix exposureCounts:
    ## E.g., replace "Signature A" with "Signature.A".
    colnames(exposureCounts) <-
      gsub(pattern = " ",replacement = ".",colnames(exposureCounts))

    ## Write attributed exposures into a SynSig formatted exposure file.
    WriteExposure(exposure_counts_dukenus,
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
