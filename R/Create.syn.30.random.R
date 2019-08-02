# This file contains functions to create "completely random"
# artificial signatures

#' Create one "random" artificial signature profile.
#'
#' @param row.names One of the \code{\link{ICAMS}} package variable such as
#'  \code{catalog.row.order[["SBS96"]]}.
#'
#' @return A single column matrix with \code{rownames} \code{row.headers} and
#'   \code{colnames} \code{"RandSig"}.
#'
#' @importFrom stats runif
#'
#' @keywords internal

CreateOneRandomMutSigProfile <- function(row.names) {
  stopifnot(!is.null(row.names))

  retval <- matrix(10^runif(length(row.names)), ncol = 1)
  # retval <- matrix(10^rnorm(length(row.names)), ncol = 1) # Too spiky
  retval <- retval / sum(retval)
  rownames(retval) <- row.names
  colnames(retval) <- "RandSig"
  return(retval)
 }

#' Create a matrix of "random" signature profiles.
#'
#' @param row.headers One of the \code{\link{ICAMS}} package variable such as
#'  \code{catalog.row.order[["SBS96"]]}.
#'
#' @param num.signatures Number of signatures to create.
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

  stopifnot(!is.null(row.headers))

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
#' @keywords internal
#'
#' @return A list with the needed synthetic means
#' and standard deviations.

CreateMeanAndStdevForSigs <-
  function(num.sigs, target.mut.mean, target.mut.sd, sig.names) {
  syn.mean <- rnorm(num.sigs, mean = target.mut.mean, sd = target.mut.sd)
  names(syn.mean) <- sig.names
  syn.sd  <- syn.mean * target.mut.sd / target.mut.mean
  names(syn.sd) <- sig.names
  return(list(syn.mean = syn.mean, syn.sd = syn.sd))
}

#' Create \code{num.exposures} signature counts
#' from a normal distribution with \code{mean} and
#' \code{sd}.
#'
#' Discard tumors with signature count = 0 or
#' signature count > \code{total.num.sigs}.
#'
#' @param num.exposures Number of exposures to create.
#'
#' @param mean Mean of distribution to draw from.
#'
#' @param sd Standard deviation of distribution to draw from.
#'
#' @param total.num.sigs Number of signatures in the "universe".
#'
#' @return A numeric vector, each element of which the number
#' of signatures in the corresponding tumor (which still remains
#' to be created.)
#'
#' @keywords internal

CreateExposuresNums <- function(num.exposures, mean,
                                sd, total.num.sigs) {
  retval <- numeric(0)
  num.exposures.to.try <- num.exposures
  while (length(retval) < num.exposures) {
    num.exposures.to.try <- 3 * num.exposures.to.try
    retval <- round(rnorm(num.exposures.to.try, mean = mean, sd = sd))
    retval <- retval[retval > 0]
    retval <- retval[retval <= total.num.sigs]
  }
  return(retval[1:num.exposures])
}

#' Select \code{num.exp} signatures from the members of \code{sig.names} and
#' create one column of an exposure matrix (as a vector).
#'
#' @param target.num.exp Number of signatures to which sample "was" exposed.
#'
#' @param all.sig.names Names of all possible signatures.
#'
#' @param target.sig.means Target means for creating mutation counts.
#'
#' @param target.sig.sds Target standard deviations for creating mutation counts.
#'
#' @return A subset of \code{sig.names} of size \code{num.exp}.
#'
#' @keywords internal

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

#' Create a pair of "random" synthetic catalogs, one for 96-channel
#' features and one for COMPOSITE features, for one set
#' of signatures.
#'
#' @param num.syn.tumors Total number of synthetic tumors to create.
#'
#' @param total.num.sigs Total number of signatures in the universe.
#'
#' @param mut.mean Mean of the log10 of the
#' number of mutations due to each signature.
#'
#' @param mut.sd Standard deviation of the log10 of
#'  the number of mutations due to each signature.
#'
#' @param num.sigs.mean Mean number of signatures contributing to each tumor.
#'
#' @param num.sigs.sd Standard deviation the number number of signatures
#' contribution to each tumor.
#'
#' @param sig.name.prefix String to put in front of an integer (as
#' string) to form an identifier for a synthetic signature.
#'
#' @param sample.name.prefix String to put in front of an integer (as
#' string) to form an identifier for a synthetic sample (tumor).
#'
#' @param composite.dir.name string indicating the name of the COMPOSITE
#' subdirectory; probably one of \code{"sa.sa.COMPOSITE"} or
#' \code{"sp.sa.COMPOSITE"}.
#'
#' @param x96.dir.name A string indicating the name of the 96-channel
#' subdirectory; probably one of \code{"sa.sa.COMPOSITE"} or
#' \code{"sp.sa.COMPOSITE"}.
#'
#' @param COMPOSITE.features Character vector containing
#' rownames for a COMPOSITE signature or catalog.
#'
#' @param overwrite If \code{TRUE} overwrite existing directories / files.
#'
#' @keywords internal

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
        ICAMS::catalog.row.order[["SBS96"]], total.num.sigs, sig.name.prefix)

    syn.COMPOSITE.sigs <-
      CreateRandomMutSigProfiles(
        COMPOSITE.features, total.num.sigs, sig.name.prefix)

    sig.info <- CreateMeanAndStdevForSigs(
      total.num.sigs, mut.mean, mut.sd, colnames(syn.96.sigs))

    buffer <- 100

    exp.nums <-
      CreateExposuresNums(
        num.exposures = num.syn.tumors + buffer,
        mean = num.sigs.mean,
        sd = num.sigs.sd,
        total.num.sigs = total.num.sigs)

    if (TRUE) {
      cat("\nCreateOneSetOfRandomCatalogs\n",
          "number of exposures per tumor statisitcs:\n")
      print(summary(exp.nums))
      cat("sd", sd(exp.nums), "\n")
    }

    exp <-
      sapply(exp.nums,
             function(x) {
               ExposureNums2Exposures(
                 x, colnames(syn.96.sigs), sig.info$syn.mean, sig.info$syn.sd) })

    test.catalog <- syn.COMPOSITE.sigs %*% exp
    stopifnot(!any(colSums(test.catalog) < 1))
    test.catalog <- round(test.catalog, digits = 0)
    zero.mutations <- colSums(test.catalog) == 0
    # colSums(test.catalog) == 0 can occur after rounding even if
    # any(colSUms(test.catalog) < 1) before rounding is FALSE, if
    # before rounding mutiple mutational classes had < 0.5 mutations.

    exp <- exp[ , !zero.mutations]
    if (ncol(exp) < num.syn.tumors)
      stop("Too many tumors with no mutations; check the code, ",
            "possibly increase the value of variable buffer")
    exp <- exp[ , 1:num.syn.tumors]
    colnames(exp) <- paste0(sample.name.prefix, 1:num.syn.tumors)

    NewCreateAndWriteCatalog(
      sigs = syn.COMPOSITE.sigs,
      exp  = exp,
      dir  = composite.dir.name,
      overwrite = overwrite)

    NewCreateAndWriteCatalog(
      sig = syn.96.sigs,
      exp = exp,
      dir = x96.dir.name,
      overwrite = overwrite)
}

#' Create a full SignatureAnalyzer / SigProfiler test data set for "random"
#' artificial signatures.
#'
#' @param top.level.dir Path to top level of directory structure to be created.
#'
#' @param num.syn.tumors Number of synthetic tumors to create.
#'
#' @param overwrite If TRUE, overwrite existing directories / files.
#'
#' @export

CreateRandomSAAndSPSynCatalogs <-
  function(top.level.dir, num.syn.tumors, overwrite = FALSE) {

  COMPOSITE.features <- c(ICAMS::catalog.row.order[["SBS1536"]],
                          ICAMS::catalog.row.order[["DBS78"]],
                          ICAMS::catalog.row.order[["ID"]])
  stopifnot(length(COMPOSITE.features) == 1697)

  if (dir.exists(top.level.dir)) {
    if (!overwrite) stop(top.level.dir, " exists and overwrite is FALSE")
  } else {
    MustCreateDir(top.level.dir)
  }

  # The following are for choosing the mean number of mutations due to each
  # synthetic signature.
  sa.mut.mean <- 2.349  # mean(log10(sa.all.real.exposures[sa.all.real.exposures >= 1]))
  sa.mut.sd   <- 0.6641 # sd(log10(sa.all.real.exposures[sa.all.real.exposures >= 1]))
  sp.mut.mean <- 2.97   # mean(log10(sp.all.real.exposures[sp.all.real.exposures >= 1]))
  sp.mut.sd   <- 0.7047 # sd(log10(sp.all.real.exposures[sp.all.real.exposures >= 1]))

  # An alternative would be:
  # per.sig.mean <- apply(sa.all.real.exposures, 1, function(x) mean(log10(x[ x > 1])
  # sa.mut.mean  <- mean(per.sig.mean, na.rm = TRUE)
  # sa.mut.sd    <- sd(per.sig.mean, na.rm = TRUE)

  sa.num.sigs.mean <- 15.525 # mean(colSums(sa.all.real.exposures > 0))
  sa.num.sigs.sd   <-  6.172 # sd(colSums(sa.all.real.exposures > 0))
  sp.num.sigs.mean <-  3.947 # mean(colSums(sp.all.real.exposures > 0))
  sp.num.sigs.sd   <-  1.331 # sd(colSums(sp.all.real.exposures > 0))

  num.sigs.to.create <- 30 # Also tired, 60 (ncol(SynSig::sa.96.sigs));
                           # this is too many.

  CreateOneSetOfRandomCatalogs(
    num.syn.tumors     = num.syn.tumors,
    total.num.sigs     = num.sigs.to.create,
    mut.mean           = sa.mut.mean,
    mut.sd             = sa.mut.sd,
    num.sigs.mean      = sa.num.sigs.mean,
    num.sigs.sd        = sa.num.sigs.sd,
    sig.name.prefix    = "SARandSig",
    sample.name.prefix = "SARandSample",
    composite.dir.name = file.path(top.level.dir, "sa.sa.COMPOSITE"), # HERE
    x96.dir.name       = file.path(top.level.dir, "sa.sa.96"),        # HERE
    COMPOSITE.features = COMPOSITE.features,
    overwrite = overwrite)

  CreateOneSetOfRandomCatalogs(
    num.syn.tumors     = num.syn.tumors,
    total.num.sigs     = num.sigs.to.create,    # was 65 -- too many! ncol(sp.sigs), # Take from actual number of SP signatures.
    mut.mean           = sp.mut.mean,
    mut.sd             = sp.mut.sd,
    num.sigs.mean      = sp.num.sigs.mean,
    num.sigs.sd        = sp.num.sigs.sd,
    sig.name.prefix    = "SPRandSig",
    sample.name.prefix = "SPRandSample",
    composite.dir.name = file.path(top.level.dir, "sp.sa.COMPOSITE"), # HERE
    x96.dir.name       = file.path(top.level.dir, "sp.sp"),           # HERE
    COMPOSITE.features = COMPOSITE.features,
    overwrite = overwrite)

  # AddAllScripts(maxK = 50)
  }

Create.syn.30.random <- function(regress = FALSE) {
  suppressWarnings(RNGkind(sample.kind="Rounding"))
  # For compatibility with R < 3.6.0
  set.seed(1443196)

  CreateRandomSAAndSPSynCatalogs("tmp.syn.30.random.sigs",
                           1000, overwrite = TRUE)
  if (regress) {
    diff.result <- Diff4SynDataSets("syn.30.random.sigs", unlink = TRUE)
    if (diff.result[1] != "ok") {
      message("\nThere was a difference, investigate\n",
              paste0(diff.result, "\n"))
    } else {
      message("\nok\n")
    }
  }
}
