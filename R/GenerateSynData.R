
#' @title Extract SynSig parameters for one mutational signature profile
#'
#' @param counts     A vector of mutation counts attributed to one signature across
#'                length(counts) samples. TODO(Steve): rename to exposures
#'
#' @param target.size The length of genomic sequence from which the counts
#'                were derived, in megabases TODO:(Steve) probably do not need this
#'
#' @return
#'   A 3-element vector with names "prevalence", "mean", and "stdev"
#'
#' @importFrom stats sd
#'
#' @export
#'
SynSigParamsOneSignature <- function(counts, target.size ) {

  prevalence <-  length(counts[counts >= 1 ]) / length(counts)

  counts.per.mb <- counts[counts >= 1 ] / target.size

  ## generate log10(mut/mb) values for mean and sd
  mean.per.mb <- mean(log10(counts.per.mb))
  sd.per.mb <- sd(log10(counts.per.mb))

  c(prob=prevalence, mean=mean.per.mb, stdev=sd.per.mb)
}


#' @title Determine 3 parameters for synthetic tumors from an exposure matrix
#'
#' @param exposures A matrix in which each column is a sample and each row is a mutation
#'         signature, with each element being the "exposure",
#'         i.e. mutation count attributed to a
#'         (sample, signature) pair.
#'
#' @param target.size Deprecated - was the length of sequence
#'         from which the counts were derived; in the future
#'         make any adjustments (e.g exome to genome) before
#'         providing the exposure matrix.
#'
#' @return A data frame with one row for
#' each of a subset of the input signatures
#'  and the following columns. Signatures not present in
#'  \code{exposures} or present only in a single tumor in
#'  \code{exposures} are removed.
#'
#' \enumerate{
#' \item the proportion of tumors with the signature
#' \item mean(log_10(mutations.per.Mb))
#' \item stdev(log_10(mutations.per.Mb))
#' }
#'
#' @export

GetSynSigParamsFromExposures <- function(exposures, target.size = 1) {

  integer.counts <- round(exposures, digits=0)
  integer.counts <- integer.counts[rowSums(integer.counts) >0 , ]
  ret1 <- apply(X=integer.counts,
                MARGIN=1,
                FUN=SynSigParamsOneSignature,
                target.size)

  # Some standard deviations can be NA (if there is only one tumor
  # with mutations for that signature). We pretend we did not see
  # these signatures. TODO(Steve): impute from similar signatures.
  if(any(is.na(ret1['stdev', ]))) {
    cat("Warning, some signatures present in only one sample, dropping:\n")
    cat(colnames(ret1)[is.na(ret1['stdev', ])], "\n")
  }
  return(ret1[,!is.na(ret1['stdev',])])
}

#' @title Write SynSig parameters --prevalence, mean(log(exposure))
#'  and sd(log(exposure)) to a file.
#'
#' @param params The parameters to write
#'
#' @param file   The path to the file to write
#'
#' @param append Whether to append to or overwrite \code{file} if it already
#'        exists.
#'
#' @importFrom utils write.table
#'
#' @export
WriteSynSigParams <- function(params, file, append = FALSE) {
  write.table(x = as.data.frame(params),
              file = file,
              sep = ",",
              col.names = ifelse(append, FALSE, NA),
              row.names = TRUE,
              append = append)
}

#' @title Create synthetic exposures based given parameters
#'
#' @return A matrix with the rows being each signature and the columns being
#' generated samples. Each entry is the count of mutations due to one
#' signature in one sample.
#'
#' @param sig.params Parameters from \code{synsig.parameters.from.attribution}
#'
#' @param num.samples Number of samples to generate
#'
#' @param name Prefix for sample identifiers in the simulated dataset
#'
#' @export

GenerateSyntheticExposures <-
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


#' @title Decide which signatures are present in the catalogs of synthetic tumors.
#'
#' @param num.tumors Number of tumors to generate
#'
#' @param prev.present Vector of prevalences, each the prevalence of 1 mutational
#'    signature
#'
#' @param sigs List(?) maybe vector(?) of signature names (?)
#'
#' @details If a tumor ends up with no signature assigned,
#' signature 1 is added as the only signature. TODO:(Steve)
#' needs redesign how to handle 0 sig assigned cases
#'
#' @keywords internal

present.sigs <-
  function(num.tumors,   ## number of tumors
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
  return(present)
}

#' @title TODO(steve) trace this to be sure we understand what it does.
#'
#' @param tumor TODO
#'
#' @param sig.interest TODO
#'
#' @param burden.per.sig TODO
#'
#' @param sd.per.sig TODO
#'
#' @details ??Determine the intensity of each
#' mutational signature in a tumor, returning mutations per mb
#' using the mean mutation burden per signature and the std dev
#'
#' @importFrom stats rbinom rnorm
#'
#' @keywords internal
get.syn.exposure <-
  function(tumor,          ## matrix with present.signatures output
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
    return(tumor)

  }

#' @title Generate synthetic spectra catalogs given
#' signature profiles and synthetic exposures.
#'
#' @param signatures The signature profiles.
#'
#' @param exposures The synthetic exposures.
#'
#' @param sample.id.suffix A string for adding a suffix to
#'  sample ID. For example, if sample.id.suffix is "abc",
#'  then SomeCancerType::s1.33 is changed to
#'  SomeCancerType::s1-abc.33. Actually, this just replaces
#'  the first "." in the sample id with "-" concatenated
#'  to sample.id.suffix. TODO(Steve): probably drop this
#'
#' @return Spectra catalog as a numeric matrix.
#'
#' @export


GenSynCatalogs <- function(signatures, exposures, sample.id.suffix = NULL) {
  exposed.sigs <- rownames(exposures)

  # It is an error if there are signatures in exposures that are not
  # in signatures.
  stopifnot(setequal(setdiff(exposed.sigs, colnames(signatures)), c()))

  # VERY IMPORTANT, the next statement guarantees that
  # the order of signatures in rows of exposures is the same as
  # the order of columns in signatures. In addition,
  # it ensure that signatures contains only signatures
  # that are present in exposures.
  #
  signatures <- signatures[ , exposed.sigs]
  catalog <- signatures %*% exposures
  i.cat <- round(catalog, digits = 0)
  if (!is.null(sample.id.suffix)) {
    newcolnames <-
      gsub(".", paste0("-", sample.id.suffix, "."),
           colnames(i.cat), fixed = TRUE)
    colnames(i.cat) <- newcolnames
  }
  return(i.cat)
  # TODO(Steve) In future, add noise
}

#' Merge 2 exposure matrices
#'
#' @param exp1 An exposure matrix
#'
#' @param exp2 An exposure matrix
#'
#' @return The column-wise merge of the two input matrices as
#' with all rownames from either matrix presevered and
#' corresponding entries filled with 0s.
#'
#' @keywords internal
Merge2Exposures <- function(exp1, exp2) {
  # Rows are signatures
  exp.m <- merge(exp1, exp2, by = 0, all = TRUE)
  rownames(exp.m) <- exp.m[ ,1]
  exp.m <- exp.m[ , -1]
  exp.m[is.na(exp.m)] <- 0
  return(as.matrix(exp.m))
}

#' Merge all exposure matrices in a list of matrices
#'
#' @param list.of.exposures A list of exposure matrices
#'
#' @return The column-wise merge of all the input matrices
#' with all rownames from all matrices presevered and
#' corresponding entries filled with 0s.
#'
#' @export

MergeExposures <- function(list.of.exposures) {
  stopifnot(length(list.of.exposures) > 0)
  if (length(list.of.exposures) == 1) {
    return(as.matrix(list.of.exposures[[1]]))
  }
  start <- list.of.exposures[[1]]
  for (i in 2:length(list.of.exposures)) {
    start <- Merge2Exposures(start, list.of.exposures[[i]])
  }
  return(as.matrix(start))
}

