# Script to run SignatureAnalyzer from the shell

library(ICAMS)  # ReadCat96, WriteCat96, etc
library(SynSig) # SignatureAnalyzerOneRun

test2 <- SignatureAnalyzerOneRun(
  signatureanalyzer.code.dir = "SignatureAnalzyer.052418/",
  input.catalog = "../vignettes/syn_3_5_40_abst_v2/sp.sp.abst.syn.data.csv",
  read.catalog.function = ReadCat96,
  out.dir = paste0("../../foo.sa.test.abst.2019.syn_3_5_40_abst_v2.02.10.work"),
  write.signature.function = WriteCat96,
  input.exposures =
    "../vignettes/syn_3_5_40_abst_v2/sp-3-5-40-abst-syn-exp.csv",
  test.only = FALSE)


