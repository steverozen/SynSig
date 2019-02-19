library(ICAMS)
library(SynSig)

set.seed(888)
INPUT <<- "INPUT_SignatureAnalyzer/"
test.abst <- SignatureAnalyzer4MatchedCatalogs(
  num.runs = 3,
  signatureanalyzer.code.dir =
    "C:/Users/steve/Documents/SynSig/data-raw/SignatureAnalzyer.052418/",
  dir.root = "c:/users/steve/Documents/0syn.3.5.40.abst.v8/",
  slice = 1,
  overwrite = TRUE,
  mc.cores = 1)

