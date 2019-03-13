# Put this file in the top level directory and start R from there.


library(SynSig)
library(ICAMS)
if (!exists("maxK.for.SA")) stop("Please set maxK.for.SA")
cat("\n\nRunning, maxK.for.SA is", maxK.for.SA, "\n\n")
RNGkind(kind = "L'Ecuyer-CMRG")
set.seed(888)

test.abst <- SignatureAnalyzer4MatchedCatalogs(
  num.runs = 20,
  signatureanalyzer.code.dir = "/home/gmssgr/bin/SignatureAnalzyer.052418/",
  dir.root = "..",
  slice = 1,
  overwrite = FALSE,
  maxK = maxK.for.SA,
  mc.cores = 20)

