# Put this file in the top level directory and start R from there.


library(SynSig)
library(ICAMS)
RNGkind(kind = "L'Ecuyer-CMRG")
set.seed(888)

test.abst <- SignatureAnalyzer4MatchedCatalogs(
  num.runs = 20,
  signatureanalyzer.code.dir = "/home/gmssgr/bin/SignatureAnalzyer.052418/",
  dir.root = "..",
  slice = 4,
  overwrite = FALSE,
  mc.cores = 20)

