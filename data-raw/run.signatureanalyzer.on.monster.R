# Put this file in the top level directory and start R from there.

library(ICAMS)
library(SynSig)

set.seed(888)
INPUT <<- "INPUT_SignatureAnalyzer/"
test.abst <- SignatureAnalyzer4MatchedCatalogs(
  num.runs = 20,
  signatureanalyzer.code.dir = "/home/gmssgr/bin/SignatureAnalzyer.052418/",
  dir.root = ".",
  slice = 3,
  overwrite = FALSE,
  mc.cores = 10)

