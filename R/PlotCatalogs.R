# PlotCatalogs.R

#' Plot the SBS96 part of a SignatureAnalyzer COMPOSITE signature or catalog
#'
#' @param catalog Catalog or signature matrix
#'
#' @param name Name of file to print to.
#'
#' @param type See \code{\link[ICAMS]{PlotCatalogToPdf}}.
#'
#' @importFrom ICAMS PlotCatSNS96ToPdf Collapse1536To96
#'
#' @keywords internal

Plot96PartOfComposite <- function(catalog, name, type = "density") {
  cat1536 <- catalog[1:1536, ]
  cat96 <- Collapse1536To96(cat1536)
  all.0 <- which(colSums(cat96) == 0)
  if (length(all.0) > 0 ) {
    cat96[ , all.0] <- 1
    cn <- colnames(cat96)
    cn[all.0] <- paste(cn[all.0], "WARNING all 0")
    colnames(cat96) <- cn
  }
  PlotCatSNS96ToPdf(catalog = cat96/sum(cat96), filename = name, type = type)
}

#' Plot the a SignatureAnalyzer COMPOSITE signature or catalog into separate pdfs
#'
#' @param catalog Catalog or signature matrix
#'
#' @param filename.header Contain path and the beginning part of the file name.
#' The name of the pdf files will be:
#' \code{filename.header}.SNS.96.pdf
#' \code{filename.header}.SNS.1536.pdf
#' \code{filename.header}.DNS.78.pdf
#' \code{filename.header}.ID.83.pdf
#'
#' @param type See \code{\link[ICAMS]{PlotCatalogToPdf}}.
#'
#' @param id A vector containing the identifiers of the samples
#' or signatures in \code{catalog}.
#'
#' @importFrom ICAMS PlotCatSNS96ToPdf Collapse1536To96
#' PlotCatSNS1536ToPdf PlotCatDNS78ToPdf PlotCatIDToPdf
#'
#' @export
PlotCatCOMPOSITE <- function(catalog, filename.header, type, id = colnames(catalog)) {

  ## Read in COMPOSITE catalogue
  test.COMPOSITE.sigs <-
    SynSig:::ReadCatCOMPOSITE(catalog)

  ## Check
  stopifnot(nrow(test.COMPOSITE.sigs) == 1697)
  # TODO WUYang: check whether the base context is in correct order

  ## Subsetting COMPOSITE catalogue
  test.SNS1536.sigs <- test.COMPOSITE.sigs[1:1536,]
  test.DNS78.sigs <- test.COMPOSITE.sigs[1537:1614,]
  test.ID83.sigs <- test.COMPOSITE.sigs[1615:1697,]

  ## Plot using ICAMS embedded plotting function


  ICAMS::PlotCatSNS1536ToPdf(test.SNS1536.sigs,
                             filename = paste0(filename.header,".SNS.1536.pdf"),
                             type = type,
                             id = id)
  ICAMS::PlotCatDNS78ToPdf(test.DNS78.sigs,
                           filename = paste0(filename.header,".DNS.78.pdf"),
                           type = type,
                           id = id)
  ICAMS::PlotCatIDToPdf(test.ID83.sigs,
                        filename = paste0(filename.header,".ID.83.pdf"),
                        type = type,
                        id = id)
}