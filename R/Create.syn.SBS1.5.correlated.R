##########################################################################################
######## Create.syn.SBS1.5.correlated.R
######## Alpha version
######## Copyright 2019 by WU Yang and Steven G ROZEN
####
#### The code is released under GPL-3
#### https://www.gnu.org/licenses/gpl-3.0.en.html
####
########
########
######## A script to generate synthetic exposures for TWO mutational signatures in
######## the log space using parameters (prob = 1, mean.log, main.stdev.log) for TWO
######## mutational signatures of interest Note: This generator generates samples
######## whose exposures of two signatures are PARTIALLY CORRELATED, This is
######## realized by letting the mean correlated.signature exposure to be the
######## exact value of main.signature exposure for each data point.
####
#### To run this code, edit parameters set starting at line 440, and then source this file
####
##########################################################################################


##########################################################################################
######## Main chunk: Function and global variables involving data generation.
##########################################################################################



##############################################################
###### Function to generate exposure for both signatures,
###### This function can force the exposure data points to be grater than
######
###### Warning: forcing dataset to have high Pearson's R^2,
###### while inputting large main.stdev.log value will cause the program
###### to reject unqualified datasets and regenerate new datasets for HOURS OR EVEN DAYS!
##############################################################

#' Function to generate exposure of two correlated signatures
#' (Example: SBS1 and SBS5) for ONE synthetic tumor.
generate.exposure.one.tumor <- function(tumor.name = "TwoCorreSigsGen::1",
                                        main.signature = "SBS5",
                                        correlated.signature = "SBS1",
                                        main.mean.log = 2.5,
                                        main.stdev.log = 0.25,
                                        correlated.stdev.log = 0.25,
                                        slope.linear = 1,
                                        main.signature.lower.thres = 50,
                                        correlated.signature.lower.thres = 30,
                                        min.main.to.correlated.ratio.linear = 1/3,  ## The lower ratio for count(SBS5) / count(SBS1) in LINEAR SPACE!
                                        max.main.to.correlated.ratio.linear = Inf)  ## The higher ratio for count(SBS5) / count(SBS1) in LINEAR SPACE!
{
  #### Matrices storing exposures for two signatures
  exposure.counts <- matrix(data = NA,              ## Referring to test.exposure.counts,
                            nrow = 2,               ## Each row refers to a signature, (Thus nrow = 2)
                            ncol = 1)               ## And each column refers to a sample. (Thus ncol = 1)


  #### Add names for signatures of interest and tumors to be generated
  rownames(exposure.counts) <- c(main.signature,correlated.signature)
  colnames(exposure.counts) <- tumor.name

  ## Repeat generating exposure for both,
  ## Until the count is greater then lower.thres AND
  ## The count(main)/count(correlated) ratio is between two thresholds
  repeat{

    ## Repeat generating the exposure for main.signature,
    ## until the count is greater than main.signature.lower.thres
    repeat{

      ## First construct a signature presence matrix for main.signature
      sig.presence <- matrix(data = c(1), nrow = 1, ncol = 1)
      rownames(sig.presence) = tumor.name
      colnames(sig.presence) = main.signature

      ## Construct mean.burden.per.sig matrix for main.signature
      mean.per.sig <- matrix(main.mean.log,1,1)
      rownames(mean.per.sig) = tumor.name
      colnames(mean.per.sig) = main.signature

      ## Construct sd.per.sig for main.signature
      sd.per.sig <- matrix(main.stdev.log,1,1)
      rownames(sd.per.sig) = tumor.name
      colnames(sd.per.sig) = main.signature

      ## Using GenerateSynExposureOneSample() to generate exposure
      ## of main.signature for a tumor.
      exposure.counts[main.signature,1] <-
        GenerateSynExposureOneSample(tumor = sig.presence,
                                     sig.interest = main.signature,
                                     burden.per.sig = mean.per.sig,
                                     sd.per.sig = sd.per.sig)
      main.signature.count.log <- log10(exposure.counts[main.signature,1])

      if(exposure.counts[main.signature,1] >= main.signature.lower.thres )
        break  ## Keep regenerating the data, until current generated exposure is greater than threshold
    }

    ## Repeat generating the exposure for correlated.signature,
    ## until the count is greater than correlated.signature.lower.thres
    repeat{

      ## Next construct a signature presence matrix for correlated.signature
      sig.presence <- matrix(data = c(1), nrow = 1, ncol = 1)
      rownames(sig.presence) = tumor.name
      colnames(sig.presence) = correlated.signature

      ## Construct mean.burden.per.sig matrix for correlated.signature
      ## To enable strong correlation between exposures of two signatures,
      ## For each tumor, the log-mean burden of correlated signature
      ## equals to log value of main signature in this tumor + log10(slope)
      mean.per.sig <- matrix(main.signature.count.log+log10(slope.linear),1,1)
      rownames(mean.per.sig) = tumor.name
      colnames(mean.per.sig) = correlated.signature

      ## Construct sd.per.sig for correlated.signature
      sd.per.sig <- matrix(correlated.stdev.log,1,1)
      rownames(sd.per.sig) = tumor.name
      colnames(sd.per.sig) = correlated.signature

      ## Using GenerateSynExposureOneSample() to generate exposure
      ## of correlated.signature for a tumor.
      exposure.counts[correlated.signature,1] <-
        GenerateSynExposureOneSample(tumor = sig.presence,
                                     sig.interest = correlated.signature,
                                     burden.per.sig = mean.per.sig,
                                     sd.per.sig = sd.per.sig)


      if( exposure.counts[correlated.signature,1] >= correlated.signature.lower.thres )
        break  ## Keep regenerating the data, until all generated exposures are greater than threshold
    }

    main.to.correlated.ratio.linear <- exposure.counts[main.signature,1] / exposure.counts[correlated.signature,1]

    if( main.to.correlated.ratio.linear >= min.main.to.correlated.ratio.linear & main.to.correlated.ratio.linear <= max.main.to.correlated.ratio.linear )
      break

  }

  return(exposure.counts)

}
