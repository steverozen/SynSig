#' Prepare the "hypermutated" segment (a.k.a "Secondary" segment
#' of a split non-hyper and hyper data set.)
#'
#' @param non.hyper.results The directory containing the the results of the analysis
#' of the non-hyper-mutated (a.k.a "PRIMARY") mutational spectra.
#'
#' @param primary.catalog The catalog of non-hyper-mutated
#' mutational spectra from which the results in \code{non.hyper.results}
#' were derived.
#'
#' @param hyper.catalog The catalog of hyper-mutated
#' mutational spectra which will be part of the input for the
#' secondary analysis.
#'
#' @param secondary.catalog The final output catalog on which
#' the secondary analysis will be performed; this is a
#' \code{cbind} of pseudo-spectra generated from the PRIMARY
#' signatures with the \code{hyper.catalog}.
#'
#' @param read.fn Function to use for reading signatures.
#'
#' @param write.fn Function to use for writing signatures.
#'
#' @param overwrite If \code{TRUE} overwrite possible previously
#' computed files and/or directories.
#'
#' @export

SignatureAnalyzerPrepHyper1Secondary <-
  function(non.hyper.results,
           primary.catalog,
           hyper.catalog,
           secondary.catalog,
           read.fn,
           write.fn,
           overwrite = TRUE) {
  if (!dir.exists(non.hyper.results)) {
    stop(non.hyper.results, "does not exist")
  }

  # Find the best run in the non-hyper-mutated data.
  best.run <-
    CopyBestSignatureAnalyzerResult(
      non.hyper.results, overwrite = overwrite)
  best.run.file <- paste0(non.hyper.results, "/best.run/sa.output.sigs.csv")
  if (!file.exists(best.run.file)) stop(best.run.file, " does not exist")
  best.sigs <- read.fn(best.run.file)

  best.exp.file <- paste0(non.hyper.results, "/best.run/sa.output.raw.exp.csv")
  if (!file.exists(best.exp.file)) stop(best.exp.file, " does not exist")
  best.exp  <- ReadExposure(best.exp.file)

  exp.sums <- rowSums(best.exp) # Sum for each signature.
  ### best.sigs not proportions, but rather are
  ### arbitrarily scaled. However, the total number of mutations
  ### in best.sigs %*% diag(exp.sums) should still approximate
  ### the orignal catalog.
  stopifnot(ncol(best.sigs) == nrow(best.exp))

  # names = TRUE handles the case where exp.sums is a scalar,
  # otherwise diag will return floor(exp.sums) X floor(exp.sums)
  # identity matrix.
  pseudo.catalog <- best.sigs %*% diag(exp.sums, names = TRUE)
  colnames(pseudo.catalog) <-
    paste0("pseudo.sample.", 1:ncol(pseudo.catalog))

  if (nrow(best.sigs) == 96) {
    # This would break on composite signatures
    PlotCatalogToPdf(
      pseudo.catalog,
      "pseudo.catalog.pdf",
      type = "counts")
  }

  # Make sure the number of mutations in pseudo.catalog is close
  # to the number in the input catalog.
  ground.truth.catalog <- read.fn(primary.catalog)
  new.mut.count <- sum(pseudo.catalog)
  orig.mut.count <- sum(ground.truth.catalog)
  cat("orig =", orig.mut.count, "new =", new.mut.count,"\n")
  pseudo.catalog <- pseudo.catalog * orig.mut.count / new.mut.count
  cat("orig =", orig.mut.count, "new =", sum(pseudo.catalog),"\n")

  hyper.catalog <- read.fn(hyper.catalog)
  hyper.catalog.plus <- cbind(hyper.catalog, pseudo.catalog)

  cat("dim of hyper.catalog.plus is", dim(hyper.catalog.plus), "\n")
  write.fn(hyper.catalog.plus, secondary.catalog)

  # We don't deal with the exposures, because we will remove the
  # pseudo-catalog from the input before assessing the extracted signatures.
  # (Must remember to do this.)

  return(NULL)

}

#' Prepare the "hypermutated" segment (a.k.a "Secondary" segment
#' of a split non-hyper and hyper data set.)
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
#' @param overwrite If \code{TRUE} overwrite existing directories and files.
#'
#' @export
#'
#' @importFrom ICAMS ReadCatalog WriteCatalog PlotCatalogToPdf

SignatureAnalyzerPrepHyper4 <-
  function(parent.dir, overwrite = FALSE) {
  if (!dir.exists(parent.dir)) stop(parent.dir, "does not exist")
  non.h.prefix <- paste0(parent.dir, "/syn.SA.hyper.low")
  h.prefix <- paste0(parent.dir, "/syn.SA.hyper.mixed")
  if (!dir.exists(non.h.prefix)) stop(non.h.prefix, "does not exist")
  if (!dir.exists(h.prefix)) stop(h.prefix, "does not exist")

  subdirs <- c("sa.sa.96", "sp.sp", "sa.sa.COMPOSITE", "sp.sa.COMPOSITE")
  read.fn <- c(ReadCatalog, ReadCatalog,
               ReadCatCOMPOSITE, ReadCatCOMPOSITE)
  write.fn <- c(WriteCatalog, WriteCatalog,
                WriteCatCOMPOSITE, WriteCatCOMPOSITE)

  # TODO(Steve): This function should be rewritten to use
  # SignatureAnalyzerPrepHyper1Secondary
  tmp.fn <- function(subdir, read.fn, write.fn) {
    dir1 <- paste0(non.h.prefix, "/", subdir)
    cat("\n\nProcessing", dir1, "\n\n")
    if (!dir.exists(dir1)) stop(dir1, "does not exist")
    dir2 <- paste0(h.prefix, "/", subdir)
    if (!dir.exists(dir2)) stop(dir2, "does not exist")
    non.hyper.results <-paste0(dir1, "/sa.results")
    if (!dir.exists(non.hyper.results)) stop(non.hyper.results, "does not exist")
    # Find the best run in the non-hyper-mutated data.
    best.run <-
      CopyBestSignatureAnalyzerResult(
        non.hyper.results, overwrite = overwrite)
    best.run.file <- paste0(dir1, "/sa.results/best.run/sa.output.sigs.csv")
    if (!file.exists(best.run.file)) stop(best.run.file, " does not exist")
    best.sigs <- read.fn(best.run.file)

    best.exp.file <- paste0(dir1, "/sa.results/best.run/sa.output.exp.csv")
    if (!file.exists(best.exp.file)) stop(best.exp.file, " does not exist")
    best.exp  <- ReadExposure(best.exp.file)

    exp.sums <- rowSums(best.exp) # Sum for each signature.
                 ### best.sigs not proportions, but rather are
                 ### arbitrarily scaled. However, the total number of mutations
                 ### in best.sigs %*% diag(exp.sums) should still approximate
                 ### the orignal catalog.
    stopifnot(ncol(best.sigs) == nrow(best.exp))

    # names = TRUE handles the case where exp.sums is a scalar,
    # otherwise diag will return floor(exp.sums) X floor(exp.sums)
    # identity matrix.
    pseudo.catalog <- best.sigs %*% diag(exp.sums, names = TRUE)
    colnames(pseudo.catalog) <-
      paste0("pseudo.sample.", 1:ncol(pseudo.catalog))

    if (nrow(best.sigs) == 96) {
      # This would break on composite signatures
      PlotCatalogToPdf(
        pseudo.catalog,
        paste0(dir2, "/pseudo.catalog.pdf"),
        type = "counts")
    }

    # Make sure the number of mutations in pseudo.catalog is close
    # to the number in the input catalog.
    ground.truth.catalog <-
      read.fn(paste0(dir1, "/ground.truth.syn.catalog.csv"))
    cat(sum(pseudo.catalog), sum(ground.truth.catalog), "\n")

    hyper.catalog.name <- paste0(dir2, "/ground.truth.syn.catalog.csv")
    hyper.catalog <- read.fn(hyper.catalog.name)
    hyper.catalog.plus <- cbind(hyper.catalog, pseudo.catalog)

    new.name <- paste0(dir2, "/prev.ground.truth.syn.catalog.csv")
    res <- file.rename(from = hyper.catalog.name, to = new.name)
    if (!res) { stop("Did not rename from ", hyper.catalog.name,
                     " to ", new.name )}

    cat("dim of hyper.catalog.plus is", dim(hyper.catalog.plus))
    write.fn(hyper.catalog.plus, hyper.catalog.name)

    # We don't deal with the exposures, because we will remove the
    # pseudo-catalog from the input before assessing the extracted signatures.
    # (Must remember to do this.)

    return(NULL)

  }

  mapply(tmp.fn, subdirs, read.fn, write.fn)

  invisible(NULL)

}

# debug(SignatureAnalyzerPrepHyper4)
# SignatureAnalyzerPrepHyper4("../syn.hyper.4.folders/")
