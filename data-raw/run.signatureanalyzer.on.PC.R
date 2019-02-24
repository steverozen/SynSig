library(ICAMS)
library(SynSig)
getwd()
set.seed(888)
INPUT <<- "INPUT_SignatureAnalyzer/"
test.abst <- SignatureAnalyzer4MatchedCatalogs(
  num.runs = 3,
  signatureanalyzer.code.dir =
    "C:/Users/rozen/Documents/SynSig/data-raw/SignatureAnalzyer.052418/",
  dir.root = "c:/users/rozen/Documents/RandomSigsV2/",
  slice = 3,
  overwrite = TRUE,
  mc.cores = 1)

