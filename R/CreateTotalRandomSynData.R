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

CreateExposuresNumsSA <- function(num.exposures) {
  mean <- 15.5252 # mean(colSums((sa.all.real.exposures > 0)))
  sd <- 6.172 # sd(colSums((sa.all.real.exposures > 0)))
  retval <- round(rnorm(3*num.exposures, mean = mean, sd = sd))
  retval <- retval[retval > 0]
  return(retval[1:num.exposures])
}

CreateExposuresNumsSP <- function(num.exposures) {
  mean <- 3.947 # mean(colSums((sp.all.real.exposures > 0)))
  sd <- 1.331 # sd(colSums((sp.all.real.exposures > 0)))
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

  num.sa.sigs <- ncol(sa.96.sigs)
  num.sp.sigs <- ncol(sp.sigs)

  syn.sa.sigs <-
    CreateRandomMutSigProfiles(
      ICAMS:::.catalog.row.order96, num.sa.sigs, "SARandSig")

  syn.sp.sigs <-
    CreateRandomMutSigProfiles(
      ICAMS:::.catalog.row.order96, num.sp.sigs, "SPRandSig")

  COMPOSITE.features <- c(ICAMS:::.catalog.row.order1536,
                          ICAMS:::.catalog.row.order.DNS.78,
                          ICAMS:::.catalog.row.order.ID)

  syn.sa.COMPOSITE.sigs <-
    CreateRandomMutSigProfiles(COMPOSITE.features, num.sa.sigs,
                               "SARandSig") # Has to be the same as for syn.sa.sigs

  syn.sp.COMPOSITE.sigs <-
    CreateRandomMutSigProfiles(COMPOSITE.features, num.sp.sigs,
                               "SPRandSig") # Has to be the same as for syn.sp.sigs

  # sp.mean <- 2.97   # mean(log10(sp.all.real.exposures[sp.all.real.exposures >= 1]))
  # sp.sd   <- 0.7047 # sd(log10(sp.all.real.exposures[sp.all.real.exposures >= 1]))
  # sa.mean <- 2.349  # mean(log10(sa.all.real.exposures[sa.all.real.exposures >= 1]))
  # sa.sd   <- 0.6641 # sd(log10(sa.all.real.exposures[sa.all.real.exposures >= 1]))

  sa.mut.mean <- 2.349  # mean(log10(sa.all.real.exposures[sa.all.real.exposures >= 1]))
  sa.mut.sd   <- 0.6641 # sd(log10(sa.all.real.exposures[sa.all.real.exposures >= 1]))
  sp.mut.mean <- 2.97   # mean(log10(sp.all.real.exposures[sp.all.real.exposures >= 1]))
  sp.mut.sd   <- 0.7047 # sd(log10(sp.all.real.exposures[sp.all.real.exposures >= 1]))

  sa.sig.info <- CreateMeanAndStdevForSigs(
    num.sa.sigs, sa.mut.mean, sa.mut.sd, colnames(syn.sa.sigs))
  sp.sig.info <- CreateMeanAndStdevForSigs(
    num.sp.sigs, sp.mut.mean, sp.mut.sd, colnames(syn.sp.sigs))

  sa.exp.nums <- CreateExposuresNumsSA(num.syn.tumors)
  names(sa.exp.nums) <- paste0("SASynSample.", 1:num.syn.tumors)
  sp.exp.nums <- CreateExposuresNumsSP(num.syn.tumors)
  names(sp.exp.nums) <- paste0("SPSynSample.", 1:num.syn.tumors)

  sa.exp <-
    sapply(sa.exp.nums,
           function(x) {
             ExposureNums2Exposures(
               x, colnames(syn.sa.sigs), sa.sig.info$syn.mean, sa.sig.info$syn.sd) })
  colnames(sa.exp) <- paste0("SASynSample.", 1:num.syn.tumors)
  sp.exp <-
     sapply(sp.exp.nums,
            function(x) {
              ExposureNums2Exposures(
                x, colnames(syn.sp.sigs), sp.sig.info$syn.mean, sp.sig.info$syn.sd) })
  colnames(sp.exp) <- paste0("SPSynSample.", 1:num.syn.tumors)
  CreateAndWriteCatalog(
    syn.sa.COMPOSITE.sigs,
    sa.exp,
    "sa.sa.COMPOSITE",
    WriteCatCOMPOSITE)

  CreateAndWriteCatalog(
    syn.sa.sigs,
    sa.exp,
    "sa.sa.96",
    WriteCat96)

  CreateAndWriteCatalog(
    syn.sp.COMPOSITE.sigs,
    sp.exp,
    "sp.sa.COMPOSITE", # Conventional name, low number of signatures per tumor,
                       # COMPOSITE features.
    WriteCatCOMPOSITE)

  CreateAndWriteCatalog(
    syn.sp.sigs,
    sp.exp,
    "sp.sp", # Conventional name low number of signatures per tumor, 96-channels.
    WriteCat96)

  }

MakeAllRandom <- function() {
  set.seed(1443196)
  CreateSAAndSPSynCatalogs("../RandomSigs", 1000, overwrite = TRUE)
}
