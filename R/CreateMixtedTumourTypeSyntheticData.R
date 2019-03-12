#' Create a full SignatureAnalyzer / SigProfiler test data set for a
#' set of different tumour type.
#'
#' @param top.level.dir Path to top level of directory structure to be created.
#'
#' @param cancer.type.strings Search the PCAWG data for tumours matching
#' these strings. Each string should identify one tumor type, for
#' some definition of tumor type. Probably the tumors in each type
#' should be non-overlapping, but the code does not enforce this and
#' does not care.
#'
#' @param num.syn.tumors Number of synthetic tumors of each cancer type
#' to create.
#'
#' @param overwrite If TRUE, overwrite existing directories / files.
#'
#' @export

CreateMixedTumorTypeSyntheticData <-
  function(top.level.dir, cancer.type.strings,
           num.syn.tumors, overwrite = FALSE) {

    SetNewOutDir(top.level.dir, overwrite)

    info.list <-
      lapply(cancer.type.strings,
             function(ca.type.str) {
               retval <-
                 SAAndSPSynDataOneCAType(
                   sa.all.real.exposures,
                   sp.all.real.exposures,
                   ca.type = ca.type.str,
                   num.syn.tumors,
                   file.prefix = "ca.type.str")
               return(retval)
             })

    # info.list has both SignatureAnalyzer and
    # SigProfiler exposures.

    sa.exposures <- lapply(info.list, function(x) x$sa.syn.exp)
    sa.exp <- MergeExposures(sa.exposures)
    sp.exposures <- lapply(info.list, function(x) x$sp.syn.exp)
    sp.exp <- MergeExposures(sp.exposures)
    # We will need the exposures later when evaluating the attributed signatures
    WriteExposure(sa.exp, OutDir("sa.exposure.csv"))
    WriteExposure(sp.exp, OutDir("sp.exposure.csv"))

    # Create synthetic mutational spectra catalogs based on SignatureAnalyzer
    # attributions

    CreateAndWriteCatalog(
      sa.COMPOSITE.sigs,
      sa.exp,
      "sa.sa.COMPOSITE",
      WriteCatCOMPOSITE)

    CreateAndWriteCatalog(
      sa.96.sigs,
      sa.exp,
      "sa.sa.96",
      WriteCatSNS96)

    # Create synthetic mutational spectra catalogs based on SigProfiler
    # attributions

    # First we need the matching between SigProfiler and SignatureAnalyzers
    # signatures.

    sp.sa.map.info <-
      MapSPToSASignatureNamesInExposure(sp.exp)

    CreateAndWriteCatalog(
      sa.COMPOSITE.sigs,
      sp.sa.map.info$exp2,
      "sp.sa.COMPOSITE",
      WriteCatCOMPOSITE)


    print(sp.sa.map.info$sp.to.sa.sig.match)

    CreateAndWriteCatalog(
      sp.sigs,
      sp.exp,
      "sp.sp",
      WriteCatSNS96)

    invisible(info.list = info.list,  sp.sa.map.info =  sp.sa.map.info)

  }

#' A simple test for CreateMixedTumorTypeSyntheticData
BladderAndUV <- function() {
  set.seed(191906)
  num.syn.tumors <- 500
  cancer.types <- c("Bladder-TCC", "Skin-Melanoma")
  retval <-
    CreateMixedTumorTypeSyntheticData(
      top.level.dir = "../Bladder-Melanoma-test",
      cancer.type.strings = cancer.types,
      num.syn.tumors = num.syn.tumors,
      overwrite = TRUE
    )
  }
