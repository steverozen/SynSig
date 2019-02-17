# Create in-memory R data from SignatureAnalyzer signature profiles.

library(ICAMS)
library(SynSig) # ReadSASigCOMPOSITE, ReadSASig96
library(usethis)

# COMPOSITE signature profiles
sa.COMPOSITE.sigs <- ReadSASigCOMPOSITE()

# 96-channel signature profiles
sa.96.sigs <-
  ReadSASig96("SignatureAnalyzer_COMPOSITE_SBS_W96.signature.031918.txt")
colnames(sa.96.sigs) <- FixSASigNames(colnames(sa.96.sigs))

usethis::use_data(sa.COMPOSITE.sigs, sa.96.sigs)
