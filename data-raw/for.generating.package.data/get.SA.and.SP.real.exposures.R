# Create in-memory R data from SignatureAnalyzer and SigProfiler attributions
# Must setwd the data-raw directory.

library(ICAMS)
library(SynSig) # FixSASigNames
library(usethis)

### Read SignatureAnalyzer exposures.

sa.no.hyper.real.exposures <-
  read.table(
    "SignatureAnalyzer_COMPOSITE_SNV.activity.no_POLE_MSI_SKIN_TMZ.031918.txt",
    sep = "\t", row.names = 1, header = TRUE)
colnames(sa.no.hyper.real.exposures) <-
  sub("__", "::", colnames(sa.no.hyper.real.exposures), fixed = TRUE)
colnames(sa.no.hyper.real.exposures) <-
  sub("_", "-", colnames(sa.no.hyper.real.exposures), fixed = TRUE)
rownames(sa.no.hyper.real.exposures) <-
  FixSASigNames(rownames(sa.no.hyper.real.exposures))
sa.no.hyper.real.exposures <- as.matrix(sa.no.hyper.real.exposures)
usethis::use_data(sa.no.hyper.real.exposures)


### Read SigProfiler exposures

sp.all.real.exposures <-
  read.csv("PCAWG_sigProfiler_SBS_signatures_in_samples.csv",
           as.is = TRUE, header = TRUE)
rownames(sp.all.real.exposures) <-
  paste0(sp.all.real.exposures$Cancer.Types, "::",
         sp.all.real.exposures$Sample.Names)
sp.all.real.exposures <- t(sp.all.real.exposures[ , -(1:3)])
sp.no.hyper.real.exposures <-
  sp.all.real.exposures[ , colnames(sa.no.hyper.real.exposures)]
sp.all.real.exposures <- as.matrix(sp.all.real.exposures)
sp.no.hyper.real.exposures <- as.matrix(sp.no.hyper.real.exposures)
usethis::use_data(sp.all.real.exposures,
                   sp.no.hyper.real.exposures)
rm(sp.all.real.exposures,
   sp.no.hyper.real.exposures,
   sa.no.hyper.real.exposures)

# Back to SA exosures, including "hypermuted" tumours.
sa.all.real.exposures <-
  read.table(
    "SignatureAnalyzer_COMPOSITE_SNV.activity.FULL_SET.031918.txt",
    sep = "\t", row.names = 1, header = TRUE)
colnames(sa.all.real.exposures) <-
  sub("__", "::", colnames(sa.all.real.exposures), fixed = TRUE)
colnames(sa.all.real.exposures) <-
  sub("_", "-", colnames(sa.all.real.exposures), fixed = TRUE)
rownames(sa.all.real.exposures) <-
  FixSASigNames(rownames(sa.all.real.exposures))
sa.all.real.exposures <- as.matrix(sa.all.real.exposures)
usethis::use_data(sa.all.real.exposures)

