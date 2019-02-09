# Script to run SignatureAnalyzer from the shell

library(SynSig)
library(data.table) # TODO(Steve) Important, Import this from ReadCat96 or update version of ICAMS in this R
options( warn = 0 )
here <- getwd()
setwd("SignatureAnalzyer.052418/") # Or the appropriate directory
INPUT <- "INPUT_SignatureAnalyzer/"
source("SignatureAnalyzer.PCAWG.function.R")
setwd(here) # This is necessary because the caller
            # as specified input and output locations
            # realtive to here.

for (i in 1:10) {
sa.ret <-
  RunSignatureAnalyzerOnFile(
    "syn_sbs3_5_40/sp.sp.RCC.and.OVA.syn.data.csv",
    ReadCat96,
    paste0("sa.test.2019.02.10.home", i),
    WriteCat96,
    test.only = FALSE)
}
