library(ICAMS)
library(SynSig)

set.seed(888)
INPUT <<- "INPUT_SignatureAnalyzer/"
test.abst <- SignatureAnalyzer4MatchedCatalogs(
  num.runs = 3,
  signatureanalyzer.code.dir = "./SignatureAnalzyer.052418/",
  dir.root = ".", slice = 1)

