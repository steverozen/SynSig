# Put this file in <top.level.dir>/sa.sa.COMPOSITE and run Rscript 3.run.SA.R
maxK.for.SA <- 30

library(SynSig)
library(ICAMS)
cat("\n\nRunning, maxK.for.SA is", maxK.for.SA, "\n\n")
RNGkind(kind = "L'Ecuyer-CMRG")
set.seed(888)

reval <- SignatureAnalyzer4MatchedCatalogs(
  num.runs = 20,
  signatureanalyzer.code.dir = "/home/gmssgr/bin/SignatureAnalzyer.052418/",
  dir.root = "..",
  slice = 3,
  overwrite = FALSE,
  maxK = maxK.for.SA,
  mc.cores = 20
  )
