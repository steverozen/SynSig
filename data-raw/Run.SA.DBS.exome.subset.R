# Script to run SignatureAnalyzer on the exome-subset of PCAWG
# DBS mutational spectra

library(ICAMS)
library(SynSig)

# num.runs                   <- 2 # 2 for debugging
num.runs                   <- 20 # 20 for production
# signatureanalyzer.code.dir <- "SignatureAnalzyer.052418" # for debugging on Laptop
signatureanalyzer.code.dir <- "/home/gmssgr/bin/SignatureAnalzyer.052418/" # for monster
input.catalog              <- "sa.DBS.exome.subset/pcawg-as-exome-DBS.csv"
read.catalog.function      <- ReadCatDNS78
out.dir                    <- "sa.DBS.exome.subset/"
write.signature.function   <- WriteCatDNS78
maxK                       <- 30
test.only                  <- FALSE
delete.tmp.files           <- TRUE
overwrite                  <- TRUE
mc.cores                   <- 20 # Will be overidden and set to 1 on Windows
verbose                    <- FALSE

# TODO(Steve): move RNGking and set.seed into SAMultiRunOneCatalog
cat("\nRunning, maxK is ", maxK, "\n\n", sep = "")
RNGkind(kind = "L'Ecuyer-CMRG")
set.seed(897)

sa.SA.DBS.exome.subset.res <-
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
