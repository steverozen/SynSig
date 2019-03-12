#' Prepare the "hypermuted" segment of a split non-hyper and hyper data set.
#'
#' @param parent.dir A directory that must contain subdirectories
#' \code{syn.SA.hyper.low} and \code{syn.SA.hyper.mixed}.
#' \code{syn.SA.hyper.low} must contain the synthetic
#' non-hypermutated data and the results of running SignatureAnalyzer on the non-hyper
#' segment, with subdirectories \code{sa.sa.96}, \code{sa.sa.COMPOSITE},
#' \code{sp.sa.COMPOSITE}, and \code{sp.sp}.
#' \code{syn.SA.hyper.mixed} must contain the synthetic hypermutated data. The results of
#' the initial SignatureAnalyzer run will be placed here to prepare this directory
#' for the second SignatureAnalyzer run.
#'
#' @param overwrite If \code{TRUE} overwrite existing directories and fies.
#'
#' @importFrom ICAMS ReadCatSNS96 WriteCatSNS96

SignatureAnalyzerPrepHyper4 <- function(parent.dir, overwrite = FALSE) {
  if (!dir.exists(parent.dir)) stop(parent.dir, "does not exist")
  non.h.prefix <- paste0(parent.dir, "/syn.SA.hyper.low")
  h.prefix <- paste0(parent.dir, "/syn.SA.hyper.mixed")
  if (!dir.exists(non.h.prefix)) stop(non.h.prefix, "does not exist")
  if (!dir.exists(h.prefix)) stop(h.prefix, "does not exist")

  subdirs <- c("sa.sa.96", "sp.sp", "sa.sa.COMPOSITE", "sp.sa.COMPOSITE")
  read.fn <- c(ReadCatSNS96, ReadCatSNS96, ReadCatCOMPOSITE, ReadCatCOMPOSITE)
  write.fn <- c(WriteCatSNS96, WriteCatSNS96, WriteCatCOMPOSITE, WriteCatCOMPOSITE)

  tmp.fn <- function(subdir, read.fn, write.fn) {
    dir1 <- paste0(non.h.prefix, "/", subdir)
    if (!dir.exists(dir1)) stop(dir1, "does not exist")
    dir2 <- paste0(h.prefix, "/", subdir)
    if (!dir.exists(dir2)) stop(dir2, "does not exist")
    non.hyper.results <-paste0(dir1, "/sa.results")
    if (!dir.exists(non.hyper.results)) stop(non.hyper.results, "does not exist")
    # Find the best run in the non-hyper-mutated data.
    best.run <- CopyBestSignatureAnalyzerResult(non.hyper.results, overwrite = overwrite)
    best.sigs <- read.fn(paste0(best.run, "/sa.output.sigs.csv"))

    best.exp  <- ReadExposure(paste0(best.run, "/sa.output.exp.csv"))
    exp.sums <- rowSums(best.exp) # Sum for each signature.
                 ### IMPORTANT, I think signatures are not proportions,
                 ### but are arbitrarily scaled, but
                 ### best.sigs %*% best.exp should still approximate
                 ### the orignal catalog.
    stopifnot(ncol(best.sigs) == nrow(best.exp))
    pseudo.catalog <- best.sigs %*% rowSums(best.exp) # This should work if
                                                      # rowSums(best.exp) is
                                                      # interpreted as a 1-column
                                                      # matrix.
    # Add check to make sure the number of mutations n pseudo.catalog is close
    # to the number in the input catalog.

    colnames(pseudo.catalog) <- paste0("pseudo.sample.", colnames(pseudo.catalog))

    hyper.catalog <- read.fn(paste0(dir2, "/ground.truth.syn.catalog"))
    hyper.catalog.plus <- cbind(hyper.catalog, pseudo.catalog)

    file.rename(from = hyper.catalog, to = paste0("prev.", hyper.catalog))
    write.fn(hyper.catalog, hyper.catalog)

    # We don't deal with the exposures, because we will remove the
    # pseudo-catalog from the input before assessing the extracted signatures.
    # (Must remember to do this.)

    return(NULL)

  }

  mapply(tmp.fn, subdirs, read.fn, write.fn)

  return(NULL)

}

# debug(SignatureAnalyzerPrepHyper4)
# SignatureAnalyzerPrepHyper4("../syn.hyper.4.folders/")
