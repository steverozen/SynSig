#' Create a full SignatureAnalyzer / SigProfiler test data set for a
#' set of various tumor types.
#'
#' @param top.level.dir Path to top level of directory structure to be created.
#'
#' @param cancer.type.strings Search the PCAWG data for tumors matching
#' these strings. Each string should identify one tumor type, for
#' some definition of tumor type. Probably the tumors in each type
#' should be non-overlapping, but the code does not enforce this and
#' does not care.
#'
#' @param num.syn.tumors Number of synthetic tumors to create
#' for each cancer type.
#'
#' @param overwrite If TRUE, overwrite existing directories / files.
#'
#' @param verbose If > 0, cat various messages.
#'
#' @export

CreateMixedTumorTypeSyntheticData <-
  function(top.level.dir, cancer.type.strings,
           num.syn.tumors, overwrite = FALSE,
           verbose = 0) {

    SetNewOutDir(top.level.dir, overwrite)

    info.list <-
      lapply(cancer.type.strings,
             function(ca.type.str) {
               if (verbose) cat("\n\nProcessing", ca.type.str, "\n\n\n")
               retval <-
                 SAAndSPSynDataOneCAType(
                   sa.all.real.exposures,
                   sp.all.real.exposures,
                   ca.type = ca.type.str,
                   num.syn.tumors,
                   file.prefix = ca.type.str)
               return(retval)
             })

    # info.list has both SignatureAnalyzer and
    # SigProfiler exposures.

    sa.exposures <- lapply(info.list, function(x) x$sa.syn.exp)
    sa.exp <- MergeExposures(sa.exposures)
    if (verbose) cat("Dimension sa.exp", dim(sa.exp), "\n")

    sp.exposures <- lapply(info.list, function(x) x$sp.syn.exp)
    sp.exp <- MergeExposures(sp.exposures)
    if (verbose) cat("Dimension sp.exp", dim(sp.exp), "\n")


    # We will need the exposures later when evaluating the attributed signatures
    WriteExposure(sa.exp, OutDir("sa.exposure.csv"))
    WriteExposure(sp.exp, OutDir("sp.exposure.csv"))

    # Create catalogs of synthetic mutational spectra
    # based on SignatureAnalyzer attributions

    CreateAndWriteCatalog(
      sa.COMPOSITE.sigs,
      sa.exp,
      "sa.sa.COMPOSITE",
      WriteCatCOMPOSITE,
      overwrite = overwrite)

    CreateAndWriteCatalog(
      sa.96.sigs,
      sa.exp,
      "sa.sa.96",
      WriteCatalog,
      overwrite = overwrite)

    # Create catalogs of synthetic mutational spectra
    # based on SigProfiler attributions

    # First we need the matching between SigProfiler and
    # SignatureAnalyzers signatures.

    sp.sa.map.info <-
      MapSPToSASignatureNamesInExposure(sp.exp)

    CreateAndWriteCatalog(
      sa.COMPOSITE.sigs,
      sp.sa.map.info$exp2,
      "sp.sa.COMPOSITE",
      WriteCatCOMPOSITE,
      overwrite = overwrite)


    if (verbose) print(sp.sa.map.info$sp.to.sa.sig.match)

    CreateAndWriteCatalog(
      sp.sigs,
      sp.exp,
      "sp.sp",
      WriteCatalog,
      overwrite = overwrite)

    AddAllScripts(maxK = 50)

    invisible(list(info.list = info.list,
                   sp.sa.map.info =  sp.sa.map.info))

  }

#' Create a specific synthetic data set of 2,700 tumors.
#'
#' @param regress If \code{TRUE}, then compare to data in \code{data-raw}
#' and report any differences; if no differences, unlink the result
#' directory.
#'
#' @export

Create.syn.many.types <- function(regress = FALSE) {
  set.seed(191906)
  num.syn.tumors <- 300 # number of tumor of each type
  top.level.dir <- "tmp.syn.many.types"

  cancer.types <- c("Bladder-TCC",    "Eso-AdenoCA",
                    "Breast-AdenoCA", "Lung-SCC",
                    "Kidney-RCC",     "Ovary-AdenoCA",
                    "Bone-Osteosarc", "Cervix-AdenoCA",
                    "Stomach-AdenoCA")
  retval <-
    CreateMixedTumorTypeSyntheticData(
      top.level.dir = top.level.dir,
      cancer.type.strings = cancer.types,
      num.syn.tumors = num.syn.tumors,
      overwrite = TRUE
    )

  if (regress) {
    if (Diff4SynDataSets("syn.many.types", unlink = TRUE) != 0) {
      cat("\nThere was a difference, investigate\n")
    } else {
      cat("\nok\n")
    }
  }

  invisible(retval)
}
