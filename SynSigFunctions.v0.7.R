## SynSigFunctions v0.7
## 2019 Jan 10
##
## Generate synthetic mutational spectra (mutations in tumours) from
##
## 1. A set of mutational signatures
##
## 2. Attribution of mutation counts in each tumour to the responsible mutational
## signatures
##
##
## Typical workflow
##
## A <- Read in matrix of attributions (signatures x samples)
## S <- Read the mutational signature profiles.
##
## P <- synsig.parameters.from.attributions (A, ....)
##
## No.noise.exposure <- target.size * generate.synthetic.exposures (P....)
##
## No.noise.specra <- ... some kind of matrix multiplicaiton ... check mSigTools/mSigAct for function
##
## Noise.spectra <- ... neg binomial resample from No.noise.spectra
##
## Use function from mSigTools to save No.noise.spectra and noise.spectra


## mSigTools for functions to read spectra and mutational signatures
# source('mSigTools.v0.11.R')

## For more information, refer to
##
## https://docs.google.com/document/d/183DzzG80BxuAdcnuL1jD1qI_YM-cgtVEdHvkbMLwHGc/edit?usp=sharing


SynSigParamsOneSignature <- function(counts, target.size ) {
    # Args:
    #   counts:      A vector of mutation counts attributed to one signature across
    #                length(counts) samples.
    #
    #   target.size: The length of genomic sequence from which the counts
    #                were dervived, in megabases
    #
    # Output:
    #   A 3-element vector with names "prevalence", "mean", and "stdv"
    #
  prevalence <-  length(counts[counts >= 1 ]) / length(counts)

  counts.per.mb <- counts[counts >= 1 ] / target.size

  ## generate log10(mut/mb) values for mean and sd
  mean.per.mb <- mean(log10(counts.per.mb))
  sd.per.mb <- sd(log10(counts.per.mb))

  c(prob=prevalence, mean=mean.per.mb, stdev=sd.per.mb)
}


synsig.params.from.attributions <- function(counts, target.size = 1) {

  ## Args:
  ##
  ## counts: A matrix in which each column is a sample and each row is a mutation
  ##         signature, with each element being the "exposure",
  ##         i.e. mutation count attributed to a
  ##         (sample, signature) pair.
  ##
  ## target.size: The length of sequence from which the counts were dervived.
  ##
  ## Output:
  ##
  ## A data frame one row for a subset of the input signatures
  ## and the following columns. Signatures not present in counts are removed.
  ##
  ## 1. the proportion of tumours with the signature
  ## 2. mean(log_10(mutations.per.Mb))
  ## 3. stddev(log_10(muations.per.Mb))

  integer.counts <- round(counts, digits=0)
  integer.counts <- integer.counts[rowSums(integer.counts) >0 , ]
  ret1 <- apply(X=integer.counts,
                MARGIN=1,
                FUN=SynSigParamsOneSignature,
                target.size)

  # Some standard deviations can be NA (if there is only one tumor
  # with mutations for that signature). We pretend we did not see
  # these signatures. TODO(steve): impute from similar signatures.
  if(any(is.na(ret1['stdev', ]))) {
    cat("Warning, some signatures present in only one sample, dropping:\n")
    cat(colnames(ret1)[is.na(ret1['stdev', ])], "\n")
  }
  return(ret1[,!is.na(ret1['stdev',])])
}

#' @title Write SynSig parameters --prevalence, mean(log(exposure))
#'  and sd(log(exposure)) to a file.
#'
#'  @param
#'    params The parameters to write
#'    file   The path to the file to write
#'    append Whether to append to or overwrite \code{file} if it already
#'        exists.
#'
WriteSynSigParams <- function(params, file, append = FALSE) {
  write.table(x = as.data.frame(params), file = file,
              sep = ",",
              col.names = NA,
              row.names = TRUE,
              append = append)
}

#' @title create synthetic exposures based given parameters
#'
#' @return A matrix with the rows being each signature and the columns being
#' generated samples. Each entry is the count of mutations due to one
#' signature in one sample.
#'
#' @param
#'   sig.params Parameters from \code{synsig.parameters.from.attribution}
#'    num.samples Number of samples to generate
#'    name Prefix for sample identifiers in the simulated dataset

generate.synthetic.exposures <-
  function(sig.params,
           num.samples = 10,
           name = 'synthetic') {

    sigs <- colnames(sig.params)
    sig.probs <- sig.params['prob', , drop = F]
    sig.burden <- sig.params['mean', , drop = F]
    sig.sd <- sig.params['stdev', , drop = F]

    sig.present <- present.sigs(num.samples,
                                sig.probs,
                                sigs)
    colnames(sig.present) <- paste(name, seq(1, num.samples), sep='.')

    ## generate a matrix of tumors with values denoting if each signature being
    ## present
    apply(sig.present,
          2,
          get.syn.exposure,
          sigs,
          sig.burden, ## burden is in mutation per megabase
          sig.sd)

  }


## present.sigs determines whether a signature is present in a tumor,
## if a tumor ends up with no signature assigned, signature 1 is added as the only signature.
### maybe needs rethinking how to handle 0 sig assigned cases
present.sigs <- function(num.tumors,   ## number of tumors
                         prev.present, ## prevalence(frequency) of a signature being present
                         sigs   ){     ## list of signatures to generate tumors with


  num.process <- length(prev.present)

  present.list <- list()

  for (i in 1:num.process){
    present.each <- rbinom(num.tumors,
                           1,
                           prob = prev.present[,i])
    present.list[[i]] <-  present.each
  }

  present <- do.call(rbind,present.list)

  rownames(present) <- sigs

  ## check if the vector has all 0s, add signature 1 to it
  for (tumor in 1:ncol(present)){
    if (all(present[,tumor] == rep(0,length(nrow(present))))){
      present['SBS1',tumor] = 1
    }
  }
  present
}

## get.syn(thetic).exposure determines the intensity of each
## mutational signature in a tumor, returning mutations per mb
## using the mean mutation burden per signature and the std dev
get.syn.exposure<- function(tumor,          ## matrix with present.signatures output
                            sig.interest,   ## signatures of interest
                            burden.per.sig, ## mutation burden in log10(muts/mb)
                            sd.per.sig      ## standard deviation of mutation burden
) {

  ## starts with individual tumors,
  active.sigs <- which(tumor != 0)

  for (sigs in active.sigs) {
    stdev <- sd.per.sig[,sigs]
    burden <- burden.per.sig[,sigs]

    ## if std dev is too big, >= 3, max = 3
    ### consider handling this different. the worry is that the variation
    ##  is too large, the sampled mutation burden will be very high,
    ### which will have a mutation burden that is not biologically possible
    if (stdev >= 3) {
      cat("Very large stdev", stdev, "\n")
      stdev = 3
    }

    ## mutational intensity follows a log normal distibution
    ## use the normal distribution with log-ed values instead
    mutations.per.mb <- 10^(rnorm(1,
                        sd = stdev,
                        mean = burden))

    tumor[sigs] <- mutations.per.mb
  }

  tumor <- as.matrix(tumor)
  names(tumor) <- sig.interest
  tumor

}
