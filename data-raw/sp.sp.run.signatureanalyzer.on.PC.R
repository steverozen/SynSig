library(ICAMS)
library(SynSig)

set.seed(888)
INPUT <<- "INPUT_SignatureAnalyzer/"
test.abst <- SignatureAnalyzer4MatchedCatalogs(
  num.runs = 3,
  signatureanalyzer.code.dir =
    "C:/Users/steve/Documents/SynSig/data-raw/SignatureAnalzyer.052418/",
  dir.root = "../mutsig-test-0.4/S.0.1.Rsq.0.1/",
  slice = 2,
  overwrite = TRUE,
  mc.cores = 1)

