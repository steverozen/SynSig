#' Run deconstructSigs attribution on a spectra catalog file.
#'
#' @param input.catalog File containing input spectra catalog. Columns are
#' samples (tumors), rows are mutation types.
#'
#' @param gt.sigs.file File containing input mutational signatures. Columns are
#' signatures, rows are mutation types.
#'
#' @param read.catalog.function Function taking a file path as
#' its only argument and returning a catalog as a numeric matrix.
#'
#' @param out.dir Directory that will be created for the output;
#' abort if it already exits.  Log files will be in
#' \code{paste0(out.dir, "/tmp")}.
#'
#' @param seed Specify the pseudo-random seed number
#' used to run deconstructSigs. Setting seed can make the
#' attribution of deconstructSigs repeatable.
#' Default: 1.
#'
#' @param maxK The maximum number of signatures to consider
#' exist in tumor spectra. Default is \code{NA}.
#' In this case the maximum would be equal to number of sigs
#' in \code{gt.sigs.file}
#'
#' @param signature.cutoff Within each tumor, signatures whose
#' relative exposures lower than this value would not be considered
#' in exposure attribution.
#' \code{signature.cutoff} was set to 0.06 by default.
#'
#' @param test.only If TRUE, only analyze the first 10 columns
#' read in from \code{input.catalog}.
#'
#' @param input.exposures A file with the synthetic exposures
#' used to generate \code{input.catalog}; if provided here,
#' this is copied over to the output directory
#' for downstream analysis.
#'
#' @param overwrite If TRUE, overwrite existing output
#'
#' @return The attributed exposure of \code{deconstructSigs}, invisibly.
#'
#' @details Creates several
#'  files in \code{paste0(out.dir, "/sa.output.rdata")}. These are
#'  TODO(Steve): list the files
#'
#' @export
#'
#' @importFrom utils capture.output


RundeconstructSigsAttributeOnly <-
  function(input.catalog,
           gt.sigs.file,
           read.catalog.function,
           out.dir,
           seed = 1,
           maxK = NA,
           signature.cutoff = 0.06,
           test.only = FALSE,
           input.exposures = NULL,
           delete.tmp.files = TRUE,
           overwrite = FALSE) {

    ## Install deconstructSigs, if not found in library.
    if("deconstructSigs" %in% rownames(installed.packages()) == FALSE){
      message("Installing deconstructSigs from CRAN...\n")

      install.packages("deconstructSigs")
    }

    ## Set seed
    set.seed(seed)
    seedInUse <- .Random.seed ## Save the seed used so that we can restore the pseudorandom series

    ## Read in spectra data from input.catalog file
    ## spectra: spectra data.frame in ICAMS format
    spectra <- read.catalog.function(input.catalog,
                                      strict = FALSE)
    if (test.only) spectra <- spectra[ , 1:10]


    ## Read in ground-truth signatures
    ## gt.sigs: signature data.frame in ICAMS format
    gt.sigs <- read.catalog.function(gt.sigs.file, strict = FALSE)

    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir)
    }

    ## Convert ICAMS-formatted spectra and signatures
    ## into deconstructSigs format
    spectra.ds <- data.frame(t(spectra))
    gt.sigs.ds <- data.frame(t(gt.sigs))

    ## Obtain attributed exposures using whichSignatures function
    ## Note: whichSignatures() can only attribute ONE tumor at each run!
    num.tumors <- nrow(spectra.ds)
    ## In each cycle, obtain attributed exposures for each tumor.
    exposures <- data.frame()

    for(ii in 1:num.tumors){
      output.list <- whichSignatures(tumor.ref = spectra.ds[ii,,drop = FALSE],
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
      exposures.one.tumor <- exposures.one.tumor * sum(spectra.ds[ii,,drop = FALSE])

      ## Bind exposures for current tumor to exposure data.frame
      exposures <- rbind(exposures,exposures.one.tumor)
    }

    ## Copy ground.truth.sigs to out.dir
    file.copy(from = gt.sigs.file,
              to = paste0(out.dir,"/ground.truth.signatures.csv"),
              overwrite = overwrite)

    ## Convert exposures to SynSig format, and output exposures.
    exposures <- t(exposures)
    write.csv(exposures,
              file = paste0(out.dir,"/attributed.exposures.csv"))

    ## Save seeds and session information
    ## for better reproducibility
    capture.output(sessionInfo(), file = paste0(dir.name,"/sessionInfo.txt")) ## Save session info
    write(x = seedInUse, file = paste0(dir.name,"/seedInUse.txt")) ## Save seed in use to a text file
    write(x = RNGInUse, file = paste0(dir.name,"/RNGInUse.txt")) ## Save seed in use to a text file

    ## Return the exposures attributed, invisibly
    invisible(exposures)
  }
