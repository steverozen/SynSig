# Script to run SignatureAnalyzer from the shell

library(ICAMS)  # Assume already installed
library(SynSig) # Assume already installed

root.dir <- "/gpfs/home/gmssgr/sa.tests/0syn.3.5.40.abst.v3/sa.sa.96/"
SignatureAnalyzerOneRun(
  signatureanalyzer.code.dir =
    "/gpfs/home/gmssgr/bin/SignatureAnalzyer.052418/",
  input.catalog = paste0(root.dir, "syn.data.csv"),
  read.catalog.function = ReadCat96,
  out.dir = paste0(root.dir, "sa.results/run.1/"),
  write.signature.function = WriteCat96)


## Assuming this script is in script.R command would be
setwd("/gpfs/home/gmssgr/sa.tests/0syn.3.5.40.abst.v3/sa.sa.96/sa.results/run.1/")
system("qsub -j -o qsub.out <pathinfo>/Rscript ./sa.sa.96.run.1.R")



