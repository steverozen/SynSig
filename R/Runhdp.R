#' Install hdp from GitHub.
#'
#' @keywords internal
Installhdp <- function(){
  message("Installing hdp from GitHub nicolaroberts/hdp ...\n")

  if ("devtools" %in% rownames(installed.packages()) == FALSE)
    install.packages("devtools")
  devtools::install_github("nicolaroberts/hdp", build_vignettes = TRUE)


}
#' Run hdp attribution on a spectra catalog file
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
#' used to run hdp. Setting seed can make the
#' attribution of hdp repeatable.
#' Default: 1.
#'
#' @param test.only If TRUE, only analyze the first 10 columns
#' read in from \code{input.catalog}.
#' Default: FALSE
#'
#' @param overwrite If TRUE, overwrite existing output.
#' Default: FALSE
#'
#' @return The attributed exposure of \code{hdp}, invisibly.
#'
#' @details Creates several
#'  files in \code{paste0(out.dir, "/sa.output.rdata")}. These are
#'  TODO(Steve): list the files
#'
#' @importFrom utils capture.output
#'
#' @export
#'
RunhdpAttributeOnly <-
  function(input.catalog,
           gt.sigs.file,
           read.catalog.function,
           out.dir,
           seedNumber = 1,
           test.only = FALSE,
           overwrite = FALSE) {
  }



#' Run hdp extraction and attribution on a spectra catalog file
#'
#' WARNING: hdp can only do exposure attribution
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
#' hdp. For a server, 30 cores would be a good
#' choice; while for a PC, you may only choose 2-4 cores.
#' By default (CPU.cores = NULL), the CPU.cores would be equal
#' to \code{(parallel::detectCores())/2}, total number of CPUs
#' divided by 2.
#'
#' @param seedNumber Specify the pseudo-random seed number
#' used to run hdp. Setting seed can make the
#' attribution of hdp repeatable.
#' Default: 1.
#'
#' @param K,K.range \code{K} is the precise value for
#' the number of signatures active in spectra (K).
#' Specify \code{K} if you know precisely how many signatures
#' are active in the \code{input.catalog}, which is the
#' \code{ICAMS}-formatted spectra file.
#'
#' \code{K.range} is A numeric vector \code{(K.min,K.max)}
#' of length 2 which tell hdp to search the best
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
#' @return The attributed exposure of \code{hdp}, invisibly.
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

Runhdp <-
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

    ## Install hdp, if not found in library
    if ("hdp" %in% rownames(installed.packages()) == FALSE)
      Installhdp()


    ## Set seed
    set.seed(seedNumber)
    seedInUse <- .Random.seed  ## Save the seed used so that we can restore the pseudorandom series
    RNGInUse <- RNGkind() ## Save the random number generator (RNG) used


    ## Read in spectra data from input.catalog file
    ## spectra: spectra data.frame in ICAMS format
    spectra <- read.catalog.function(input.catalog,
                                     strict = FALSE)
    if (test.only) spectra <- spectra[ , 1:10]
    convSpectra <- t(spectra)

    number.channels <- dim(spectra)[1]
    number.samples <- dim(spectra)[2]

    ## Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
    } else {
      dir.create(out.dir, recursive = T)
    }

    ## Determine the best number of signatures (K.best).
    ##
    ## hdp accepts an initial guess of number of signatures (K.initial), and later
    ## determine the best number of signatures (K.best)
    ##
    ## If K is provided, use K as the K.initial.
    ## If K.range is provided, use the largest value as the K.initial.
    ##
    ##
    if(bool1)
      K.initial <- K
    if(bool2)
      K.initial <- max(K.range)

    print(paste0("Assuming there are ",K.initial," signatures active in input spectra."))
    print(paste0("But the final number of signatures may not equal to ",K.initial,"\n."))

    ## Run hdp main program.
    ## Step 1: initialize hdp object
    {
      ## initialise hdp
      ppindex <- c(0, 1, rep(2,number.samples))
      cpindex <- c(1, 2, rep(3,number.samples))

      hdpObject <- hdp::hdp_init(ppindex = ppindex,
                      cpindex = cpindex,
                      hh = rep(1,number.channels),
                      alphaa = rep(1,3), alphab = rep(1,3))
      num.process <- hdp::numdp(hdpObject)

      ## Add data
      hdpObject <- hdp::hdp_setdata(hdpObject, 3:num.process, convSpectra)

      hdp::numdp(hdpObject)

      ## hdp also has to enter number of signatures in advance, but the final result doesn't necessarily equal to the initial input value
      ## When the number of sample is too large, hdp may eat up all your RAMs!
      hdpObject <- hdp::dp_activate(hdpObject,
                                    1:num.process,
                                    K.initial,
                                    seed=seedNumber)	## K.initial initial components(start with 30 signatures)

      hdpObject

      ## Release the occupied RAM by dp_activate
      gc()
      gc()
      gc()
    }

    ## Step 2: run 4 independent sampling chains
    {
      ## Run four independent posterior sampling chains
      chlist <- vector("list", 4)	#4 is too much here!

      i <- 1 ## Should execute manually rather than for cycle, because it is too slow!
      for (i in 1:4){
        chlist[[i]] <- hdp::hdp_posterior(hdpObject,
                                     burnin=4000,	## 4000 is too large, but necessary
                                     n=50,	## n here refers to the times of posterior sampling after burnin. To be faster, n can be set to 50.
                                     space=50,	## space is the time of iterations between two samplings. In this case, I need to iterate 9000 times.
                                     cpiter=3,
                                     seed=seedNumber)
      }

      ## Generate the original multi_chain for the sample
      mut_example_multi <- hdp::hdp_multi_chain(chlist)
      mut_example_multi
    }

    ## Step 3: Plot the diagnostics of sampling chains.
    {
      ## Plotting using hdp functions
      ## Plotting device on the server does not work
      ## Need to plot the file into a pdf



      ## Draw the DP oscillation plot for mut_example_multi(original_sample)
      {
        pdf(file = paste0(out.dir,"/original_sample.pdf"))

        par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
        p1 <- lapply(hdp::chains(mut_example_multi),
                     hdp::plot_lik, bty="L", start=500)
        p2 <- lapply(hdp::chains(mut_example_multi),
                     hdp::plot_numcluster, bty="L")
        p3 <- lapply(hdp::chains(mut_example_multi),
                     hdp::plot_data_assigned, bty="L")

        dev.off()
      }

      ## Extract components(here is signature) with cosine.merge = 0.90 (default)
      mut_example_multi_extracted <- hdp::hdp_extract_components(mut_example_multi)
      mut_example_multi_extracted


      ## Generate a pdf for mut_example_multi_extracted
      {
        pdf(file = paste0(out.dir,"/signature_hdp_embedded_func.pdf"))
        ## Draw the DP oscillation plot for mut_example_multi_extracted
        par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
        p1 <- lapply(hdp::chains(mut_example_multi_extracted),
                     hdp::plot_lik, bty="L", start=500)
        p2 <- lapply(hdp::chains(mut_example_multi_extracted),
                     hdp::plot_numcluster, bty="L")
        p3 <- lapply(hdp::chains(mut_example_multi_extracted),
                     hdp::plot_data_assigned, bty="L")

        ## Draw the computation size plot
        par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
        hdp::plot_comp_size(mut_example_multi_extracted, bty="L")

        ## Close the PDF device so that the plots are exported to PDF
        dev.off()
      }
    }

    ## Step 4: Using hdp samples to extract signatures
      {

        ## Calculate the mutation composition in each signature:
        extractedSignatures <- hdp::comp_categ_distn(mut_example_multi_extracted)$mean
        dim(extractedSignatures)
        ## Ans is 10(Components) 96(Mutational types) for hdp_500_10
        ## Ans is 9(Components) for hdp_500_12
        extractedSignatures <- t(extractedSignatures)

        ## Output the signatures extracted
        write.catalog.function(x = extractedSignatures,
                               file = paste0(out.dir,"/extracted.signatures.csv"))
      }


    ## Step 5: Using hdp samples to attribute exposure counts.
    {

      ## Calculate the exposure probability of each signature(component) for each tumor sample(posterior sample corresponding to a dirichlet process node):
      ## This is the probability distribution of signatures(components) for all tumor samples(DP nodes)
      ## exposureProbs proves to be the normalized signature exposure all 100 tumor samples

      exposureProbs <- hdp::comp_dp_distn(mut_example_multi_extracted)$mean
      dim(exposureProbs)
      exposureProbs <- exposureProbs[3:dim(exposureProbs)[1],]
      rownames(exposureProbs) <- rownames(convSpectra)[1:dim(exposureProbs)[1]]
      colnames(exposureProbs) <- colnames(extractedSignatures)[5:dim(extractedSignatures)[2]]
      dim(exposureProbs)


      ## Calculate signature exposure counts from signature exposure probability
      ## Unnormalized exposure counts = Normalized exposure probability * Total mutation count in a sample
      sample_mutation_count <- apply(convSpectra,1,sum)

      exposureCounts <- matrix(nrow = dim(exposureProbs)[1], ncol = dim(exposureProbs)[2])
      dimnames(exposureCounts) <- dimnames(exposureProbs)
      for (sample in seq(1,dim(exposureProbs)[1])) exposureCounts[sample,] <- sample_mutation_count[[sample]] * exposureProbs[sample,]

      ## Reshape exposureCounts matrix for better compatibility with plot_signatures.R
      exposureCounts <- cbind(rownames(exposureCounts),exposureCounts)
      colnames(exposureCounts)[1] <- "Samples"

      ## Next, write the exposureCounts to a file
      ## Write attributed exposures into a SynSig formatted exposure file.
      WriteExposures(exposureCounts,
                     paste0(out.dir,"/attributed_exposures.csv"))
    }



      ## Save seeds and session information
      ## for better reproducibility
      capture.output(sessionInfo(), file = paste0(out.dir,"/sessionInfo.txt")) ## Save session info
      write(x = seedInUse, file = paste0(out.dir,"/seedInUse.txt")) ## Save seed in use to a text file
      write(x = RNGInUse, file = paste0(out.dir,"/RNGInUse.txt")) ## Save seed in use to a text file

      ## Return a list of signatures and exposures
      invisible(list("signature" = extractedSignatures,
                     "exposure" = exposureCounts))
  }
