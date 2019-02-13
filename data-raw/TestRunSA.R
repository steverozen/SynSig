# Script to run SignatureAnalyzer from the shell

wrapper <- function() {
  library(ICAMS) # ReadCat96, WriteCat96, etc
  library(SynSig)
  # library(data.table) # TODO(Steve) Important, Import this from ReadCat96 or update version of ICAMS in this R
  options( warn = 0 )
  here <- getwd()
  setwd("SignatureAnalzyer.052418/") # Or the appropriate directory
  INPUT <- "INPUT_SignatureAnalyzer/"
  source("SignatureAnalyzer.PCAWG.function.R")
  setwd(here) # This is necessary because the caller
  # as specified input and output locations
  # realtive to here.

  for (i in 1:5) {
      RunSignatureAnalyzerOnFile(
        # "syn_sbs3_5_40/sp.sp.RCC.and.OVA.syn.data.csv",
        "../vignettes/syn_3_5_40_abst_v2/sp.sp.abst.syn.data.csv",
        ReadCat96,
        paste0("../../saX.test.abst.2019.syn_3_5_40_abst_v2.02.10.work", i),
        # paste0("../../sa.test.2019.02.10.home", i),
        WriteCat96,
        input.exposures =
          "../vignettes/syn_3_5_40_abst_v2/sp-3-5-40-abst-syn-exp.csv",
        test.only = FALSE)
  }


}

test2 <- SignatureAnalyzerOneRun(
  signatureanalyzer.code.dir = "SignatureAnalzyer.052418/",
  input.catalog = "../vignettes/syn_3_5_40_abst_v2/sp.sp.abst.syn.data.csv",
  read.catalog.function = ReadCat96,
  out.dir = paste0("../../foo.sa.test.abst.2019.syn_3_5_40_abst_v2.02.10.work"),
  write.signature.function = WriteCat96,
  input.exposures =
    "../vignettes/syn_3_5_40_abst_v2/sp-3-5-40-abst-syn-exp.csv",
  test.only = FALSE)

# Run this from data-raw
debug(RunSignatureAnalyzerOnFile)
test3 <- SignatureAnalyzer4MatchedCatalogs(
  num.runs = 20,
  signatureanalyzer.code.dir = "SignatureAnalzyer.052418/",
  dir.root = "../../0syn.3.5.40.abst.v3/",
  test.only = TRUE)

