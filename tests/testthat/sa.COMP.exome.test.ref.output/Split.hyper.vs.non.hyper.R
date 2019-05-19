# Split the SignatureAnalyzer exome-subset COMPOSITE spectra into separate mutational
# spectra catalogs for hyper-mutated and non-hyper-mutated spectra based on Jaegil's
# previous classification.

library(SynSig)

setwd("C:/Users/steve/Documents/SynSig/data-raw/sa.COMPOSITE.exome.subset/")

non.hyper.ids <- colnames(sa.no.hyper.real.exposures)

hyper.ids <-
  setdiff(colnames(sa.all.real.exposures), non.hyper.ids)

all.spectra <- ReadCatCOMPOSITE("pcawg-as-exome-COMPOSITE.csv")
new.colnames <- sub("..SP", "::SP", colnames(all.spectra), fixed = TRUE)
new.colnames2 <- sub(".", "-", new.colnames, fixed = TRUE)
colnames(all.spectra) <- new.colnames2

hyper.spectra <- all.spectra[ , hyper.ids]

non.hyper.spectra <- all.spectra[ , non.hyper.ids]

stopifnot(ncol(hyper.spectra) + ncol(non.hyper.spectra) == ncol(all.spectra))

WriteCatCOMPOSITE(hyper.spectra, "hyper-pcawg-as-exome-COMPOSITE.csv")

WriteCatCOMPOSITE(non.hyper.spectra, "non-hyper-pcawg-as-exome-COMPOSITE.csv")
