# Create in-memory R data from SigProfiler signature 96 changel signature profile.

library(ICAMS)
library(SynSig) # ReadSASigCOMPOSITE, ReadSASig96
library(usethis)

sp.sigs <-
  ReadCat96("sigProfiler_SBS_signatures_2018_03_28.csv",
            strict = FALSE)

usethis::use_data(sp.sigs)
