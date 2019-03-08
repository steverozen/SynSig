CatalogMutationBurdenHist <- function(path) {
  cat <- ReadCatSNS96(path)
  invisible(
    hist(log10(colSums(cat))

         ))
}

