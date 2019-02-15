# Test scripts to run SignatureAnalyzer

library(ICAMS)
library(SynSig)


test3 <- SignatureAnalyzer4MatchedCatalogs(
  num.runs = 20,
  signatureanalyzer.code.dir = "SignatureAnalzyer.052418/",
  dir.root = "../../0syn.3.5.40.abst.v3/", # top level of data / results tree
  slice = 2) # slice 2 is sp.sp

