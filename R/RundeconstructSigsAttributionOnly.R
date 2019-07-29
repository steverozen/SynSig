#' Install deconstructSigs from CRAN
InstalldeconstructSigs <- function(){
  message("Installing deconstructSigs from CRAN...\n")
  install.packages("deconstructSigs")
}

#' Run deconstructSigs attribution on a spectra catalog file
#' and known signatures.
#'
#' @param input.catalog File containing input spectra catalog. Columns are
#' samples (tumors), rows are mutation types.
#'
#' @param gt.sigs.file File containing input mutational signatures. Columns are
#' signatures, rows are mutation types.
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
#' used to run deconstructSigs. Setting seed can make the
#' attribution of deconstructSigs repeatable.
#' Default: 1.
#'
#' @param test.only If TRUE, only analyze the first 10 columns
#' read in from \code{input.catalog}.
#' Default: FALSE
#'
#' @param overwrite If TRUE, overwrite existing output.
#' Default: FALSE
#'
#' @return The attributed exposure of \code{deconstructSigs}, invisibly.
#'
#' @details Creates several
#'  files in \code{paste0(out.dir, "/sa.output.rdata")}. These are
#'  TODO(Steve): list the files
#'
#' @importFrom utils capture.output
#'
#' @export

RundeconstructSigsAttributeOnly <-
  function(input.catalog,
           gt.sigs.file,
           read.catalog.function,
           out.dir,
           seedNumber = 1,
           test.only = FALSE,
           overwrite = FALSE) {

    ## Install deconstructSigs, if not found in library.
    if("deconstructSigs" %in% rownames(installed.packages()) == FALSE)
      InstalldeconstructSigs()


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
    ## gt.sigs: signature data.frame in ICAMS format
    gtSignatures <- read.catalog.function(gt.sigs.file)

    ## Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }

    ## Convert ICAMS-formatted spectra and signatures
    ## into deconstructSigs format
    convSpectra <- data.frame(t(spectra))
    gt.sigs.ds <- data.frame(t(gtSignatures))

    ## Obtain attributed exposures using whichSignatures function
    ## Note: deconstructSigs::whichSignatures() can only attribute ONE tumor at each run!
    num.tumors <- nrow(convSpectra)
    ## In each cycle, obtain attributed exposures for each tumor.
    exposures <- data.frame()

    for(ii in 1:num.tumors){
      output.list <- deconstructSigs::whichSignatures(tumor.ref = convSpectra[ii,,drop = FALSE],
                                                      signatures.ref = gt.sigs.ds,
                                                      contexts.needed = TRUE)
      ## names(output.list): [1] "weights" "tumor"   "product" "diff"    "unknown"
      ## $weights: attributed signature exposure (in relative percentage)
      ## Note: sum of all exposure may be smaller than 1
      ## $tumor: input tumor spectrum
      ## $product: Reconstructed catalog = product of signatures and exposures
      ## = $weights %*% gt.sigs.ds
      ## $diff: $product - $tumor
      ## $unknown: 100% - $weights
      ## (percentage of exposures not attributed by this program)

      ## Obtain absolute exposures for current tumor
      exposures.one.tumor <- output.list$weights
      exposures.one.tumor <- exposures.one.tumor * sum(convSpectra[ii,,drop = FALSE])

      ## Bind exposures for current tumor to exposure data.frame
      exposures <- rbind(exposures,exposures.one.tumor)
    }



    ## Write exposure counts in ICAMS and SynSig format.
    exposures <- t(exposures)
    WriteExposure(exposures,
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
    invisible(exposures)
  }
