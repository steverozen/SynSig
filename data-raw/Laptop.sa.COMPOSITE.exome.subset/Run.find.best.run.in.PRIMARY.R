# Find the best run in sa.COMPOSITE.exome.subset/PRIMARY/
# Run this from data-raw/sa.COMPOSITE.exome.subset/.
# The folder sa.COMPOSITE.exome.subset/PRIMARY/ should
# have subdirectories sa.run.1, ..., sa.run.<n>

library(SynSig)

chosen.dir <-
  CopyBestSignatureAnalyzerResult(
    sa.results.dir = "PRIMARY/",
    overwrite = TRUE,
    verbose = TRUE
  )
