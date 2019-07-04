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


#' Generate correlated exposures for one tumor
#'
#' Function to generate exposure of two correlated signatures
#' (Example: SBS1 and SBS5) for ONE synthetic tumor.
#'
#' Warning 1: forcing dataset to have high Pearson's R^2,
#' while inputting large main.stdev.log value will cause
#' the program to reject unqualified datasets and regenerate
#' new datasets for HOURS OR EVEN DAYS!
#'
#' Warning 2: For replication of paperresults, only specify
#' parameters as indicated by supplementary information!
#'
#' @param tumor.name Name of synthetic tumor you want to generate.
#' (default: "TwoCorreSigsGen::1")
#'
#' @param main.signature Name of a signature with smaller variance
#' in the log10 space. (default: "SBS5")
#'
#' @param correlated.signature Name of a signature with larger variance
#' in the log10 space. (default: "SBS1")
#'
#' @param main.mean.log Mean of log10(mutation burden of \code{main.signature})
#'
#' @param main.stdev.log Standard deviation of log10(mutation burden
#' of \code{main.signature})
#'
#' @param correlated.stdev.log Contribute to part of the standard deviation of
#' log10(mutation burden of correlated.signature). In this script, the s.d. of
#' log10(mutation burden of correlated.signature)
#' = main.stdev.log + correlated.stdev.log
#'
#' @param slope.linear Average ratio of mutation burden of \code{correlated.signature}
#' over mutation burden of \code{main.signature}
#'
#' @param main.signature.lower.thres Minimum mutation burden (number of mutations)
#' induced by \code{main.signature} in each tumor.
#'
#' @param correlated.signature.lower.thres Minimum mutation burden (number of mutations)
#' induced by \code{correlated.signature} in each tumor.
#'
#' @param min.main.to.correlated.ratio.linear Minimum ratio of
#' \code{main.signature} over mutation burden of
#' \code{correlated.signature} in each tumor.
#'
#' @param max.main.to.correlated.ratio.linear Maximum ratio of
#' \code{main.signature} over mutation burden of
#' \code{correlated.signature} in each tumor.
#'
#' @keywords internal
#'
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


## Generic function for generating signature exposure matrix.
## sample.number specify the number of tumors in this matrix
## Parameters control tumors' exposure amount to the main.signature(SBS5) and correlated signature(SBS1)
## This function is defined for tuning the parameters, and reshape the datasets more flexibily
## Signature exposure matrix can be served as "ground-truth" or "benchmark" attribution,
## and thus can be used to compare with attribution results from different software tools.

#' Generate correlated exposures for multiple tumors
#'
#' Wrapper function around generate.exposure.one.tumor():
#' A function to generate exposure of two correlated signatures
#' (Example: SBS1 and SBS5) for \code{sample.number} (e.g. 500) synthetic tumors.
#'
#' NOTE: \code{pearson.r.2.lower.thres} and \code{pearson.r2.higher.thres}
#' are used to constraint the Pearson's R^2 of mutation burdens of two signatures
#' in multiple tumors.
#'
#'
#' @param main.signature Name of a signature with smaller variance
#' in the log10 space. (Default: "SBS5")
#'
#' @param correlated.signature Name of a signature with larger variance
#' in the log10 space. (Default: "SBS1")
#'
#' @param sample.number Number of tumors whose mutation burdens
#' will be generated. (Default: 500)
#'
#' @param name.prefix Prefix of tumor name.
#' (Default: \code{"TwoCorreSigsGen"})
#' By default, the name of tumors to be created will be:
#' TwoCorreSigGen::1, TwoCorreSigGen::2, TwoCorreSigGen::3...
#'
#' @param main.mean.log Mean of log10(mutation burden of \code{main.signature})
#'
#' @param main.stdev.log Standard deviation of log10(mutation burden
#' of \code{main.signature})
#'
#' @param correlated.stdev.log Contribute to part of the standard deviation of
#' log10(mutation burden of correlated.signature). In this script, the s.d. of
#' log10(mutation burden of correlated.signature)
#' = main.stdev.log + correlated.stdev.log
#'
#' @param slope.linear Average ratio of mutation burden of \code{correlated.signature}
#' over mutation burden of \code{main.signature}
#'
#' @param main.signature.lower.thres Minimum mutation burden
#' (number of mutations) induced by \code{main.signature} in each tumor.
#'
#' @param correlated.signature.lower.thres Minimum mutation burden
#' (number of mutations) induced by \code{correlated.signature} in each tumor.
#'
#' @param pearson.r.2.lower.thres Minimum Pearson's R^2 of mutation burdens
#' of two signatures in \code{sample.number} tumors.
#'
#' @param pearson.r.2.higher.thres Maximum Pearson's R^2 of mutation burdens
#' of two signatures in \code{sample.number} tumors.
#'
#' @param min.main.to.correlated.ratio.linear Minimum ratio of
#' \code{main.signature} over mutation burden of
#' \code{correlated.signature} in each tumor.
#'
#' @param max.main.to.correlated.ratio.linear Maximum ratio of
#' \code{main.signature} over mutation burden of
#' \code{correlated.signature} in each tumor.
#'
#'
#' @export
#'
generate.exposure <- function(main.signature,
                              correlated.signature,
                              sample.number,
                              name.prefix = "TwoCorreSigsGen",
                              main.mean.log = 2.5,
                              main.stdev.log = 0.25,
                              correlated.stdev.log = 0.25,
                              slope.linear = 1,
                              main.signature.lower.thres = 50,
                              correlated.signature.lower.thres = 30,
                              pearson.r.2.lower.thres = 0.1,
                              pearson.r.2.higher.thres = 1.0,
                              min.main.to.correlated.ratio.linear = 1/3,    ## The lower ratio for count(SBS5) / count(SBS1) in LINEAR SPACE!, minimum value is 0
                              max.main.to.correlated.ratio.linear = Inf) ## The higher ratio for count(SBS5) / count(SBS1) in LINEAR SPACE! maximum value is Inf
{
  #### Matrices storing exposures for two signatures
  #### For a specific tumor sample, exposure.mb records exposure number of a signature per unit megabase
  #### exposure.counts records exposure number of a signature throughout whole genome/exome
  #### We use exposure.counts to test the performance of differnet tools.
  exposure.counts <- matrix(data = NA,              ## Referring to test.exposure.counts,
                            nrow = 2,               ## Each row refers to a signature, (Thus nrow = 2)
                            ncol = sample.number)   ## And each column refers to a sample. (Thus ncol = sample.number)


  #### Add names for signatures of interest and tumors to be generated
  rownames(exposure.counts) <- c(main.signature,correlated.signature)
  colnames(exposure.counts) <- paste(name.prefix,seq(1,sample.number),sep = "::") ## Designate names for synthetic tumors

  repeat{
    #### Repeat generating exposures for both signatures
    #### Note: For each individual tumor, we need to discard the generated exposure,
    #### if the count(main.signature) < main.signature.lower.thres
    for( ii in 1:sample.number ){ ## Generate a single tumor at a time in order to let the discarding process faster
      temp.exposure.counts <-
        generate.exposure.one.tumor(tumor.name = paste(name.prefix,"::",ii),
                                    main.signature = main.signature,
                                    correlated.signature = correlated.signature,
                                    main.mean.log = main.mean.log,
                                    main.stdev.log = main.stdev.log,
                                    correlated.stdev.log = correlated.stdev.log,
                                    slope.linear = slope.linear,
                                    main.signature.lower.thres = main.signature.lower.thres,
                                    correlated.signature.lower.thres = correlated.signature.lower.thres,
                                    min.main.to.correlated.ratio.linear = min.main.to.correlated.ratio.linear,   ## The lower ratio for count(SBS5) / count(SBS1) in LINEAR SPACE!
                                    max.main.to.correlated.ratio.linear = max.main.to.correlated.ratio.linear)

      exposure.counts[,ii] <- temp.exposure.counts
    }

    a <- exposure.counts[main.signature,]
    b <- exposure.counts[correlated.signature,]
    ## If Pearson's R^2 > pearson.r.2.lower.thres in log-transformed dataset, then accept this dataset
    if( pearson.r.2.higher.thres >= cor(log10(a),log10(b)) ^ 2  & cor(log10(a),log10(b)) ^ 2 > pearson.r.2.lower.thres )
      break
  }

  return(exposure.counts)
}
##############################################################
###### Function to plot the scatterplot between exposures of two signatures
##############################################################
plot.correlation.scatterplot <-
  function (x,y,
            xlab = NULL,
            ylab = NULL,
            main = NULL,
            optional.remarks = "",
            ...)
{
  ## Check whether x and y are vectors of same length
  stopifnot(is.vector(x) & is.vector(y))
  stopifnot(length(x) == length(y))

  plot(x = x, y = y,
       xlab = xlab,
       ylab = ylab,
       main = main,
       ... = ...)

  mtext(paste("Pearson R^2 = ", round(cor(x,y)^2,3), "; ",
              "x.stdev = ", round(sqrt(var(x)),3), "; ",
              "y.stdev = ", round(sqrt(var(y)),3), "; ",
              sep = ""),
        cex = 0.8) ## By default, the mtext() function adds text at the top margin
  mtext(paste("x.mean = ", round(mean(x),3), "; ",
              "y.mean = ", round(mean(y),3), "; ",
              "Number of data points: ", length(x),
              sep = ""),
        line = 1,
        cex = 0.8) ## By default, the mtext() function adds text at the top margin
  mtext(optional.remarks, ## Adds a line of optional.remarks above the previous text line.
        line = 2,
        cex = 0.8) ## By default, the mtext() function adds text at the top margin
}

plot.correlation.scatterplot.from.objects <-
  function(pdf.filename,
           parameter.obj,
           exposure.counts,
           xlim=c(0,4),
           ylim=c(0,4),
           ...)
{

  main.signature <- parameter.obj$main.signature
  correlated.signature <- parameter.obj$correlated.signature
  slope.linear <- parameter.obj$slope.linear

  pdf(pdf.filename)

  plot.correlation.scatterplot(x = log(exposure.counts[main.signature,], base = 10),
                               y = log(exposure.counts[correlated.signature,], base = 10),
                               xlab = paste("log10( ",main.signature," )",sep = ""),
                               ylab = paste("log10( ",correlated.signature," )",sep = ""),
                               main = "",
                               optional.remarks = paste("slope.linear = ", slope.linear,"; ",
                                                        sep = ""),
                               xlim = xlim,
                               ylim = ylim,
                               ... = ...)
  dev.off()
}
##########################################################################################
######## Main chunk: Generate Dataset
######## the exposure of two signatures are correlated.
######## Generate the exposure for SBS5 (main.signature) first,
######## Then generate SBS1 (correlated.signature) exposure correlated to the SBS5 exposure.
##########################################################################################
######## If you want to customize the dataset's Pearson R^2,
######## you need to change the standard deviations of two signatures.
######## i.e., main.stdev.log and correlated.stdev.log.
##########################################################################################

#' Example function for generating SBS1-SBS5-correlated
#' Synthetic data.
CreateSBS1SBS5CorrelatedSyntheticData <-
  function(){

    #### Load SigProfiler signature profiles
    {
      sp.omics.sigs <- list()   ## A list storing SigProfiler genome/exome/flat signatures
      sp.signature.types <- c("SBS.96","SBS.1536","DBS.78","ID.83")


      sp.omics.sigs[["SBS.96"]] <-  get.signatures(signature.file = './SP.signatures.attributions/sigProfiler_SBS_signatures.csv', ## This is a human genome SBS signature
                                                   exome.op = .h19.96.sureselect.v6.op)  ## This is an opportunity stored in mSigTools, for exome (Because the sum of total occurence is 57Mb)
      ## exome.op is required by get.signatures, but we do not use the exome signature -- only the genome signature.

    }

    #### Load SignatureAnalyzer signatures; not used at present
    if (FALSE) {
      sa.omics.sigs <- list()
      sa.signature.types <- c("COMPOSITE.SBS.96","COMPOSITE.SBS.1536","DBS.78","ID.83")

      sa.omics.sigs[["COMPOSITE.SBS.96"]] <-  get.signatures(signature.file = './SA.signatures.attributions/SignatureAnalyzer_COMPOSITE_SBS_W96.signature.031918.txt', ## This is a human genome SBS signature
                                                             exome.op = .h19.96.sureselect.v6.op)
    }




    #### To Run the code, edit parameters from here down and then source this file
    {

      dataset.name <- "S.0.5.Rsq.0.3" ## The dataset.name encodes the parameters for the synthetic data, but this is just a convention
      ## S.0.5 refers to slope.linear = 0.5
      ## Rsq.0.3 refers to control the Pearson's R^2 from 0.29 to 0.31
      ## i.e., pearson.r.2.lower.thres = 0.29 and pearson.r.2.higher.thres = 0.31
      cat(paste("Specifying dataset.name as: ",dataset.name,"...\n",sep = ""))
      ## make a directory to store the dataset.
      dir.name <- paste("./Synthetic_datasets/", dataset.name,sep = "")

      ## We do not want to overwite an existing data set
      if(dir.exists(dir.name)){
        stop(paste("Folder ",dir.name," already exists!\n",sep = ""))
      }

      dir.create(dir.name, recursive = TRUE) ## Create ./Synthetic_datasets/ folder if current working directory doesn't have one
      cat(paste("Output folder for this dataset is: ",dir.name,"\n",sep = ""))

      dataset <- list() # This will contain the data set and the parameters used to generate it

      ## WARNING:
      ## Exposure generation function will repeat generating exposure counts using mean and stdev parameters,
      ## until the dataset has a Pearson's R^2 which falls between two boundaries of Pearson's R^2.
      ## Below are a group of parameters which have been tested successfully.
      ## If you intend to lower the Pearson's R^2, do remember to increase the main.stdev.log and correlated.stdev.log.
      ## Otherwise, the exposure generation will keep generating and discarding datasets!
      dataset$parameter <-
        list(main.signature = "SBS5",           ## The name of the main signature whose exposure can vary freely.
             correlated.signature = "SBS1",     ## The name of the correlated signature whose exposure is influenced by and co-varies with the exposure of main.signature.
             ## In this study, it defaults as "SBS1"
             sample.number = 500,       ## The number of synthetic tumors you want to generate
             main.mean.log = 2.5,       ## The mean of log(count(SBS5),base = 10)
             main.stdev.log = 0.3,      ## The standard deviation of log(count(SBS5),base = 10)
             correlated.stdev.log = 0.4,    ## The ADDED standard deviation of log(count(SBS1),base = 10).
             ## This parameter is ADDED stdev because
             ## based on the mechanism to generate the count, log10(count(SBS1)) inherently has a stdev = slope * main.stdev.log
             slope.linear = 0.5,    ## The ratio for: (Correlated exposure) / (Main exposure) IN LINEAR SPACE!
             main.signature.lower.thres = 100,      ## This program will force the exposure count of
             ## main.signature to be greater than this threhold.
             correlated.signature.lower.thres = 1,      ## This program will force the exposure count of
             ## correlated.signature to be greater than this threhold.
             pearson.r.2.lower.thres = 0.29,      ## Lower boundary of Pearson's R^2
             pearson.r.2.higher.thres = 0.31,     ## Upper boundary of Pearson's R^2
             min.main.to.correlated.ratio.linear = 1/3, ## The lower ratio for count(SBS5) / count(SBS1) in LINEAR SPACE!
             max.main.to.correlated.ratio.linear = Inf) ## The higher ratio for count(SBS5) / count(SBS1) in LINEAR SPACE!)

      #### Generate exposure matrices for main.signature
      cat("Generating ground-truth exposures according to parameters specified...\n")
      dataset$exposure <- generate.exposure.from.parameter.obj(dataset$parameter)

      #### Plot out the scatter plot for the two correlated exposures
      cat("Plotting correlation scatterplot for exposures of two signatures...\n")
      plot.correlation.scatterplot.from.objects(pdf.filename = paste(dir.name,"/Scatterplot.",dataset.name,".pdf",sep = ""),
                                                parameter.obj = dataset$parameter,
                                                exposure.obj = dataset$exposure,
                                                xlim = c(0,4),
                                                ylim = c(0,4))


      #### Using the exposure count generate and plot synthetic catalog with or without noise.
      dataset$spectra <- generate.spectra(exposure.object = dataset$exposure,
                                          signature.matrix = sp.omics.sigs[["SBS.96"]][["exome"]][,c(dataset$parameter$main.signature,dataset$parameter$correlated.signature)])

      plot.spectra(spectra.object = dataset$spectra,
                   filename = paste(dir.name,"/Spectra.plot.",dataset.name,".pdf",sep = ""))
      cat("Spectra generated.")

      #### Output Duke-NUS formatted exposure catalogs and exposure.counts

      spectra.to.txt(spectra.object = dataset$spectra,filename = paste(dir.name,"/",dataset.name,sep = ""))
      exposure.to.txt(exposure.object = dataset$exposure,filename = paste(dir.name,"/",dataset.name,sep = ""))
      spectra.to.csv(spectra.object = dataset$spectra,filename = paste(dir.name,"/",dataset.name,sep = ""))
      exposure.to.csv(exposure.object = dataset$exposure,filename = paste(dir.name,"/",dataset.name,sep = ""))

      #### Output parameters used for better reproducibility
      write.table(t(data.frame(dataset$parameter)),paste(dir.name,"/",dataset.name,".parameters.txt",sep = ""),col.names = F,quote = F,sep = ":\t")

      #### Save the session image
      save.image(paste(dir.name,"/",dataset.name,".RData",sep = ""))

      cat(paste("Result files stored in folder ",dir.name,"\n",sep = ""))
      cat("\n\nData generation has been finished successfully!\n\n")
    }

  }
