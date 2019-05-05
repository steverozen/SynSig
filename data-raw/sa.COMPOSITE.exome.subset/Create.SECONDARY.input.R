# Prepare the SECONDARY subdirecgtory of data-raw/sa.COMPOSITE.exome.subset/.
# Run this from data-raw/sa.COMPOSITE.exome.subset/.
# The folder sa.COMPOSITE.exome.subset/PRIMARY/ should
# have subdirectories sa.run.1, ..., sa.run.<n>

library(SynSig)

SignatureAnalyzerPrepHyper1Secondary(
  non.hyper.results = "PRIMARY",
  primary.catalog   = "non-hyper-pcawg-as-exome-COMPOSITE.csv",
  hyper.catalog     = "hyper-pcawg-as-exome-COMPOSITE.csv",
  secondary.catalog = "SECONDARY.catalog.csv",
  read.fn           = ReadCatCOMPOSITE,
  write.fn          = WriteCatCOMPOSITE)

