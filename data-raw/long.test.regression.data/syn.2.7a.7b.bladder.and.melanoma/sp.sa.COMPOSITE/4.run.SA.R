# Put this file in <top.level.dir>/sp.sa.COMPOSITE and run Rscript 4.run.SA.R
maxK.for.SA <- 30

library(SynSig)
library(ICAMS)
cat("

Running, maxK.for.SA is", maxK.for.SA, "

")
RNGkind(kind = "L'Ecuyer-CMRG")
set.seed(888)

reval <- SignatureAnalyzer4MatchedCatalogs(
  num.runs = 20,
  signatureanalyzer.code.dir = "/home/gmssgr/bin/SignatureAnalzyer.052418/",
  dir.root = "..",
  slice = 4,
  overwrite = FALSE,
  maxK = maxK.for.SA,
  mc.cores = 20
  )
