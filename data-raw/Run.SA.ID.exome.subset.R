# This script runs against an implicit argument

library(ICAMS)
library(SynSig)

num.runs                   <- 20 # 2 for debugging
# signatureanalyzer.code.dir <- "SignatureAnalzyer.052418" # for debugging on Laptop
signatureanalyzer.code.dir <- "/home/gmssgr/bin/SignatureAnalzyer.052418/"
input.catalog              <- "sa.ID.exome.subset/pcawg-as-exome-ID.csv"
read.catalog.function      <- ReadCatID
out.dir                    <- "sa.ID.exome.subset/"
write.signature.function   <- WriteCatID
maxK                       <- 30
test.only                  <- FALSE
delete.tmp.files           <- TRUE
overwrite                  <- TRUE
mc.cores                   <- 20 # Will be overidden and set to 1 on Windows
verbose                    <- TRUE

# TODO(Steve): move RNGking and set.seed into SAMultiRunOneCatalog
cat("\nRunning, maxK is ", maxK, "\n\n", sep = "")
RNGkind(kind = "L'Ecuyer-CMRG")
set.seed(897)

sa.SA.ID.exome.subset.res <-
  SAMultiRunOneCatalog(
    num.runs                   = num.runs,
    signatureanalyzer.code.dir = signatureanalyzer.code.dir,
    input.catalog              = input.catalog,
    read.catalog.function      = read.catalog.function,
    out.dir                    = out.dir,
    write.signature.function   = write.signature.function,
    maxK                       = 30,
    tol                        = 1e-7,
    test.only                  = test.only,
    delete.tmp.files           = delete.tmp.files,
    overwrite                  = overwrite,
    mc.cores                   = mc.cores,
    verbose                    = verbose)
