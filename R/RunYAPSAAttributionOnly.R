#' Run YAPSA attribution on a spectra catalog file.
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
#' @param input.exposures A file with the synthetic exposures used to generate
#' \code{input.catalog}; if provided here,
#' this is copied over to the output directory
#' for downstream analysis.
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
           maxK = NA,
           signature.cutoff = 0.06,
           test.only = FALSE,
           input.exposures = NULL,
           delete.tmp.files = TRUE,
           overwrite = FALSE) {

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

    ###### Exposure attribution: given gt.sigs (ground truth signatures) 
    ###### and  (tumor spectra),
    ###### determine the exposures of these signatures in tumor spectra

    in_signatures_df <- gt.sigs  ## Known signature matrix

    #### Tumor spectra matrix and related parameters
    in_mutation_catalogue_df <- spectra ## Converted spectra matrix
    size <- colSums(in_mutation_catalogue_df) ## Total mutation count of each spectrum

    ## Plotting parameter - maximum height in the plot
    ymax <- rep(0.4,ncol(in_mutation_catalogue_df))
    names(ymax) <- colnames(in_mutation_catalogue_df)

    #### Using Linear Combination Decomposition to attribute exposures
    #### YAPSA::LCD() is not recommended. The author recommended YAPSA::LCD_complex_cutoff()
	if(FALSE){
    LCD_object <- YAPSA::LCD(in_mutation_catalogue_df, in_signatures_df,
                      in_per_sample_cutoff = 0)  ## Signature presence cutoff is set to 0
      ## This means that as long as SBS1 and SBS5 are attributed to >0% of total mutations,
      ## they will be kept as candidate signatures.
      ## This makes sense because in this synthetic dataset,
      ## all tumors have exposures to SBS1 and SBS5.

      ## YAPSA::LCD() returns a data.frame with attributed exposures.
      ## Note: For each tumor spectrum, sum of attributed exposures by YAPSA::LCD()
      ## does not equal to the sum of ground-truth exposures!
      class(LCD_object) ## [1] "data.frame"
      dim(LCD_object) ## [1] 2 500
    }
	
    ## YAPSA also supports different presence cutoff for different signatures,
    ## this is done by providing different values of cutoff in LCD_complex_cutoff function.
    ## Authors suggest to use YAPSA::LCD_complex_cutoff() rather than YAPSA::LCD() in most cases.
    LCD_complex_object <- YAPSA::LCD_complex_cutoff(in_mutation_catalogue_df,
                                             in_signatures_df,
                                             in_cutoff_vector = c(0,0), ## If there are 2 signatures in the spectra,
                                             ## you must provide a
                                             in_rescale = TRUE)  ## Rescale signature exposures so that the sum of exposure for each tumor
    ## equals to the exposure sum in original spectra
    ## This prevents the difference between original spectra and observed spectra
    class(LCD_complex_object) ## [1] "list"
    names(LCD_complex_object) ## For detail, see YAPSA user manual
    ##[1] "exposures"                   "norm_exposures"
    ##[3] "signatures"                  "choice"
    ##[5] "order"                       "residual_catalogue"
    ##[7] "rss"                         "cosDist_fit_orig_per_matrix"
    ##[9] "cosDist_fit_orig_per_col"    "sum_ind"
    ##[11] "out_sig_ind_df"              "aggregate_exposures_list"

    ## Exposures generated by LCD_complex_object()
    ## does not equal to exposures generated by LCD()
    dim(LCD_complex_object$exposures) == dim(LCD_object) ## [1] TRUE
    LCD_complex_object$exposures == LCD_object ## [1] FALSE

    ## For each tumor spectrum, $exposures (the exposure counts attributed by LCD_complex_object())
    ## sums up to the total mutation counts in 500 tumors in the dataset.
    ## But $norm_exposures (relative exposure probs attributed by LCD_complex_object())
    ## sums up to 500 only.
    sum(LCD_complex_object$exposures) == sum(spectra) ## [1] TRUE
    sum(LCD_complex_object$norm_exposures) ## [1] 500

    ## For each tumor spectrum, sum of normalized attributed exposures by LCD_complex_cutoff()
    ## does not equal to the sum of ground-truth exposures.
    all( colSums(LCD_complex_object$norm_exposures) == colSums(spectra) ) ## [1] FALSE

    ## Export attributed exposure probs
    LCD_exposure_prob <- LCD_complex_object$norm_exposures
    ## Export attributed exposure counts
    LCD_exposure_counts <- LCD_complex_object$exposures ## Export exposure probs


    ## Output exposure counts
    LCD_exposure_counts_output <- data.frame("Samples" = colnames(LCD_exposure_counts),t(LCD_exposure_counts))
    write.table(x = LCD_exposure_counts_output,file = paste(out.dir,"attributed_exposures.txt",sep="/"),sep = "\t", row.names = FALSE, quote = FALSE)


    ######################
    ###### Save workspace, command history,
    ###### seeds and session information for better reproducibility
    ######################
    save.image(paste(out.dir,"/","YAPSA_",dataset.name,".RData",sep=""))  ## Save workspace image to an .RData file
    capture.output(sessionInfo(), file = paste(out.dir,"other_output_from_YAPSA","sessionInfo.txt",sep="/")) ## Save session info
    write(x = seedInUse, file = paste(out.dir,"other_output_from_YAPSA","seedInUse.txt",sep="/")) ## Save seed in use to a text file
    savehistory(file =  paste(out.dir,"/","YAPSA_",dataset.name,".Rhistory",sep="")) ## Save history file to an .Rhistory file

  }
