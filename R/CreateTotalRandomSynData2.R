# Create "completely random" signatures


#' Create one "random" signature profile.
#'
#' @param row.headers One of the \code{\link{ICAMS}} package variable such as
#'  \code{\link{.catalog.row.order96}}.
#'
#' @return A single column matrix with rownames \code{row.headers} and
#'   colname \code{RandSig}.
#'
#' @keywords internal

CreateOneRandomMutSigProfile <- function(row.headers) {
  retval <- matrix(10^runif(length(row.headers)), ncol = 1)
  # retval <- matrix(10^rnorm(length(row.headers)), ncol = 1) # Too spiky
  retval <- retval / sum(retval)
  rownames(retval) <- row.headers
  colnames(retval) <- "RandSig"
  return(retval)
 }

#' Create a matrix of "random" signature profiles.
#'
#' @param row.headers One of the \code{\link{ICAMS}} package variable such as
#'  \code{\link{.catalog.row.order96}}.
#'
#' @param num.signatures Number of signaures to create.
#'
#' @param sig.name.prefix The signatures will be named \code{<sig.name.prefix>1},
#' \code{<sig.name.prefix>2}, etc.
#'
#' @return A \code{num.signatures}-column
#'   matrix with rownames \code{row.headers}.
#'
#' @keywords internal

CreateRandomMutSigProfiles <-
  function(row.headers, num.signatures, sig.name.prefix) {
  retval <- lapply(1:num.signatures,
                function (x) CreateOneRandomMutSigProfile(row.headers))
  retval <- as.matrix(data.frame(retval))
  colnames(retval) <- paste0(sig.name.prefix, 1:num.signatures)
  return(retval)
}

#' Create means and standard deviations of log10 mutation counts for synthetic
#' signatures.
#'
#' @param num.sigs Number of signatures for which means and standard deviations
#' are needed.
#'
#' @param target.mut.mean Target number of mean of log10 of number
#'  of mutations per tumor due to one signature.
#'
#' @param target.mut.sd Target of standard deviation of  the the log10 of the
#' number of mutations per tumor due to one signature.
#'
#' @param sig.names Names for the output vectors.
#'
#' @return A list with the needed synthetic means and standard deviations.

CreateMeanAndStdevForSigs <-
  function(num.sigs, target.mut.mean, target.mut.sd, sig.names) {
  syn.mean <- rnorm(num.sigs, mean = target.mut.mean, sd = target.mut.sd)
  names(syn.mean) <- sig.names
  syn.sd  <- syn.mean * target.mut.sd / target.mut.mean
  names(syn.sd) <- sig.names
  return(list(syn.mean = syn.mean, syn.sd = syn.sd))
}

CreateExposuresNums <- function(num.exposures, mean, sd) {
  retval <- round(rnorm(3*num.exposures, mean = mean, sd = sd))
  retval <- retval[retval > 0]
  return(retval[1:num.exposures])
}

#' Select \code{num.exp} from the elemtns of \code{sig.names}.
#'
#' @param target.num.exp Number of signatures to which sample "was" exposed.
#'
#' @param all.sig.names Names of all possible signatures.
#'
#' @param target.sig.means Target means for creating mutation counts.
#'
#' @param target.sig.sds Target standard deviations for creating mutation counts.
#'
#' @retun A subset of \code{sig.names} of size \code{num.exp}.
#'
ExposureNums2Exposures <-
  function(target.num.exp, all.sig.names, target.sig.means, target.sig.sds) {
    stopifnot(names(target.sig.means) == names(target.sig.sds))
    stopifnot(all.sig.names == names(target.sig.means))
    all.len <- length(all.sig.names)
    index.to.use <- sample.int(all.len, size=target.num.exp, replace = FALSE)
    retval <- numeric(all.len) # all zeros
    retval[index.to.use] <-
      10^rnorm(target.num.exp,
               mean = target.sig.means[index.to.use],
               sd = target.sig.sds[index.to.use])
    names(retval) <- all.sig.names
    return(retval)
  }

CreateOneSetOfRandomCatalogs <-
  function(num.syn.tumors,
           total.num.sigs,
           mut.mean,
           mut.sd,
           num.sigs.mean,
           num.sigs.sd,
           sig.name.prefix,
           sample.name.prefix,
           composite.dir.name,
           x96.dir.name,
           COMPOSITE.features,
           overwrite = FALSE) {

    syn.96.sigs <-
      CreateRandomMutSigProfiles(
        ICAMS:::.catalog.row.order96, total.num.sigs, sig.name.prefix)

    syn.COMPOSITE.sigs <-
      CreateRandomMutSigProfiles(
        COMPOSITE.features, total.num.sigs, sig.name.prefix)

    sig.info <- CreateMeanAndStdevForSigs(
      total.num.sigs, mut.mean, mut.sd, colnames(syn.96.sigs))

    exp.nums <-
      CreateExposuresNums(num.syn.tumors, num.sigs.mean, num.sigs.sd)

    exp <-
      sapply(exp.nums,
             function(x) {
               ExposureNums2Exposures(
                 x, colnames(syn.96.sigs), sig.info$syn.mean, sig.info$syn.sd) })
    colnames(exp) <- paste0(sample.name.prefix, 1:num.syn.tumors)

    CreateAndWriteCatalog(
      syn.COMPOSITE.sigs,
      exp,
      composite.dir.name,
      WriteCatCOMPOSITE,
      overwrite = overwrite)

    CreateAndWriteCatalog(
      syn.96.sigs,
      exp,
      x96.dir.name,
      WriteCat96,
      overwrite = overwrite)
}

#' Create a full SignatureAnalyzer / SigProfiler test data set.
#'
#' @param top.level.dir Path to top level of directory structure to be created.
#'
#' @param num.sig.tumors Number of synthetic tumors to create.
#'
#' @param overwrite If TRUE, overwrite existing directories / files.
#'
#' @export

CreateSAAndSPSynCatalogs <-
  function(top.level.dir, num.syn.tumors, overwrite = FALSE) {
  SetNewOutDir(top.level.dir, overwrite)

  COMPOSITE.features <- c(ICAMS:::.catalog.row.order1536,
                          ICAMS:::.catalog.row.order.DNS.78,
                          ICAMS:::.catalog.row.order.ID)
  sa.mut.mean <- 2.349  # mean(log10(sa.all.real.exposures[sa.all.real.exposures >= 1]))
  sa.mut.sd   <- 0.6641 # sd(log10(sa.all.real.exposures[sa.all.real.exposures >= 1]))
  sp.mut.mean <- 2.97   # mean(log10(sp.all.real.exposures[sp.all.real.exposures >= 1]))
  sp.mut.sd   <- 0.7047 # sd(log10(sp.all.real.exposures[sp.all.real.exposures >= 1]))

  sa.num.sigs.mean <- 15.525 # mean(colSums((sa.all.real.exposures > 0)))
  sa.num.sigs.sd   <-  6.172 # sd(colSums((sa.all.real.exposures > 0)))
  sp.num.sigs.mean <-  3.947 # mean(colSums((sp.all.real.exposures > 0)))
  sp.num.sigs.sd   <-  1.331 # sd(colSums((sp.all.real.exposures > 0)))


  CreateOneSetOfRandomCatalogs(
    num.syn.tumors     = num.syn.tumors,
    total.num.sigs     = ncol(sa.96.sigs), # Take from actual number of SA signatures.
    mut.mean           = sa.mut.mean,
    mut.sd             = sa.mut.sd,
    num.sigs.mean      = sa.num.sigs.mean,
    num.sigs.sd        = sa.num.sigs.sd,
    sig.name.prefix    = "SARandSig",
    sample.name.prefix = "SARandSample",
    composite.dir.name = "sa.sa.COMPOSITE",
    x96.dir.name       = "sa.sa.96",
    COMPOSITE.features = COMPOSITE.features,
    overwrite = overwrite)

  CreateOneSetOfRandomCatalogs(
    num.syn.tumors     = num.syn.tumors,
    total.num.sigs     = ncol(sp.sigs), # Take from actual number of SP signatures.
    mut.mean           = sp.mut.mean,
    mut.sd             = sp.mut.sd,
    num.sigs.mean      = sp.num.sigs.mean,
    num.sigs.sd        = sp.num.sigs.sd,
    sig.name.prefix    = "SPRandSig",
    sample.name.prefix = "SPRandSample",
    composite.dir.name = "sp.sa.COMPOSITE",
    x96.dir.name       = "sp.sp",
    COMPOSITE.features = COMPOSITE.features,
    overwrite = overwrite)
  }

MakeAllRandom <- function() {
  set.seed(1443196)
  CreateSAAndSPSynCatalogs("C:/Users/steve/Documents/RandomSigsV2/",
                           1000, overwrite = TRUE)
}
