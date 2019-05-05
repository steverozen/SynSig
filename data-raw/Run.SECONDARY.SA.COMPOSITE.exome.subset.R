# Script to run SignatureAnalyzer on the exome-subset of PCAWG
# COMPOSITE mutational spectra, *non-hyper-mutated specta only*.
# Jaegil calls this "PRIMARY". Run this from SynSig/data-raw/

library(ICAMS)
library(SynSig)

laptop.test <- FALSE

if (laptop.test) {
  input.catalog              <- "Laptop.sa.COMPOSITE.exome.subset/SECONDARY.catalog.csv"
  out.dir                    <- "Laptop.sa.COMPOSITE.exome.subset/SECONDARY/"
  test.only                  <- TRUE
  num.runs                   <- 2
  maxK                       <- 10
  signatureanalyzer.code.dir <- "C:/Users/steve/Documents/SynSig/data-raw/SignatureAnalzyer.052418/"
  stopifnot(getwd() == "C:/Users/steve/Documents/SynSig/data-raw")
} else {
  input.catalog              <- "sa.COMPOSITE.exome.subset/SECONDARYcatalog.csv"
  out.dir                    <- "sa.COMPOSITE.exome.subset/SECONDARY/"
  test.only                  <- FALSE
  num.runs                   <- 20
  maxK                       <- 80
  signatureanalyzer.code.dir <- "/home/gmssgr/bin/SignatureAnalzyer.052418/"
}

read.catalog.function      <- SynSig:::ReadCatCOMPOSITE
write.signature.function   <- WriteCatCOMPOSITE
delete.tmp.files           <- TRUE
overwrite                  <- TRUE
mc.cores                   <- 20 # Will be overidden and set to 1 on Windows
verbose                    <- FALSE

# TODO(Steve): move RNGking and set.seed into SAMultiRunOneCatalog
cat("\nRunning, maxK is ", maxK, "\n\n", sep = "")
RNGkind(kind = "L'Ecuyer-CMRG")
set.seed(892)

sa.COMPOSITE.exome.subset.res <-
  SAMultiRunOneCatalog(
    num.runs                   = num.runs,
    signatureanalyzer.code.dir = signatureanalyzer.code.dir,
    input.catalog              = input.catalog,
    read.catalog.function      = read.catalog.function,
    out.dir                    = out.dir,
    write.signature.function   = write.signature.function,
    maxK                       = maxK,
    tol                        = 1e-7,
    test.only                  = test.only,
    delete.tmp.files           = delete.tmp.files,
    overwrite                  = overwrite,
    mc.cores                   = mc.cores,
    verbose                    = verbose)
