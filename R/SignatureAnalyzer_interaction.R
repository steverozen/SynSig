# SignatureAnalyzer_interaction.R

library(data.table)

#####################################################################
# Functions to read SignatureAnalyzer signature catalogs from disk
# into standard in-memory representations.
#####################################################################

#' Functions to canonicalize mutation type names
#'
#' \code{CanonicalizeSAx} For 96 and 1536 SNSs.
#'
#' \code{CanonicalizeSADBS} For 78 DNS (doublet substitutions).
#'
#' \code{CanonicalizeSAID} For indels.
#'
#' @return Mutation type names in canonical in-memory format.
#'
#' @name CanonicalizeSAMutTypes
NULL

#' @rdname CanonicalizeSAMutTypes
#' @keywords internal
CanonicalizeSAx <- function(feature.names) {
  xx <-
    as.matrix(data.frame(strsplit(feature.names,
                                  "_at_",
                                  fixed = TRUE)))
  ref.gt.var       <- xx[1, ]
  before.ref.after <- xx[2, ]
  var <- substring(ref.gt.var, 3, 3)
  return(paste0(before.ref.after, var))
}

#' @rdname CanonicalizeSAMutTypes
#' @keywords internal
CanonicalizeSADBS <- function(feature.names) {
  return(sub(">", "", feature.names, fixed = TRUE))
}

#' @rdname CanonicalizeSAMutTypes
#' @keywords internal
CanonicalizeSAID <- function(feature.names) {
  return(gsub("_", ":", feature.names, fixed = TRUE))
}

#' Read a catalog of spectra or signatures in SignatureAnalyzer format.
#'
#' @param path File to read
#'
#' @param canonicalize.rows A function to get the rownames and put the rows
#'                      in the correct order.
#'
#' @param exp.nrow          Expected number of rows. Stop if the actual
#'                      number of rows in path is different.
#'
#' @param row.order         Canonical row order for these mutation types.
#'
#' @return   A catalog of spectra or signatures in canonical in-memory format.
#'
#' @importFrom data.table fread
#'
#' @keywords internal
ReadAnySACat <- function(path, canonicalize.rows, exp.nrow, row.order) {
  dt <- as.data.frame(fread(path))
  stopifnot(nrow(dt) == exp.nrow)
  feature <- which(names(dt) == "feature")
  feature.names <- dt[ , feature]
  dt <- dt[ , -1 * feature]

  rownames(dt) <- canonicalize.rows(feature.names)

  out <- dt[row.order, ]

  return(as.matrix(out))
}

#' Read Catalog Functions
#'
#' Read a catalog or signature in SignatureAnalyzer format from path
#'
#' \code{ReadSACat96} Read a 96 SNS catalog from path
#'
#' \code{ReadSACat1536} Read a 1536 SNS catalog from path
#'
#' \code{ReadSACatDBS} Read a 78 DNS catalog from path
#'
#' \code{ReadSACatID} Read a ID (insertion/deletion) catalog from path
#' Please take note that deletion Repeat Size ranges from 0 to 5+
#' in the catalog, but for plotting and end user documentation
#' it ranges from 1 to 6+.
#'
#' @param path The file to read from.
#'
#' @return A catalog in canonical in-memory format.
#'
#' @name ReadSACat
NULL

#' @rdname ReadSACat
#' @export
ReadSACat96 <- function(path) {
  ReadAnySACat(path, CanonicalizeSAx, 96, ICAMS:::.catalog.row.order96)
}

#' @rdname ReadSACat
#' @export
ReadSACat1536 <- function(path) {
  ReadAnySACat(path, CanonicalizeSAx, 1536,
               ICAMS:::.catalog.row.order1536)
}

#' @rdname ReadSACat
#' @export
ReadSACatDBS <- function(path) {
  ReadAnySACat(path, CanonicalizeSADBS, 78,
               ICAMS:::.catalog.row.order.DNS.78)
}

#' @rdname ReadSACat
#' @export
ReadSACatID <- function(path) {
  ReadAnySACat(path, CanonicalizeSAID, 83,
               ICAMS:::.catalog.row.order.ID)
}

#' Read SignatureAnalyzer COMPOSITE signatures from fixed files
#'
#' @return A matrix (data.frame?) of COMPOSITE signatures
#'
#' @export
ReadSASigsCOMPOSITE <- function() {

  sa1536.file <-
    "SignatureAnalyzer_COMPOSITE_SBS_W1536.signature.031918.txt"
  saDBS.file  <-
    "SignatureAnalyzer_COMPOSITE_DBS.signature.042018.txt"
  saID.file   <-
    "SignatureAnalyzer_COMPOSITE_ID.signature.031918.txt"

  sa1536.sig  <- ReadSACat1536(sa1536.file)
  saDBS.sig   <- ReadSACatDBS(saDBS.file)
  saID.sig    <- ReadSACatID(saID.file)

  colnames(sa1536.sig) <-
    sub("_SNV", "", colnames(sa1536.sig), fixed = TRUE)
  colnames(saDBS.sig)  <-
    sub("_DNP", "", colnames(saDBS.sig), fixed = TRUE)
  colnames(saID.sig)   <-
    sub("_INDEL", "", colnames(saID.sig), fixed = TRUE)
  stopifnot(colnames(sa1536.sig) == colnames(saDBS.sig))
  stopifnot(colnames(sa1536.sig) == colnames(saID.sig))
  return(rbind(sa1536.sig, saDBS.sig, saID.sig))
}

#####################################################################
# Convenience function, used for dumping COMPOSITE catalogs
# to disk without any fancy re-writing of the representation of the
# mutational signatures.
#####################################################################

#' Probably not needed
#'
#' @param df A data.frame
#'
#' @param file File to write to
#'
#' @param rowname.name Name of first column in file
#'
#' @importFrom data.table fwrite as.data.table

fwriteDataFrame <- function(df, file, rowname.name = "mutation.type") {
  df <- as.data.frame(df)
  df[ , rowname.name] <- rownames(df)
  df <- df[ , c(ncol(df), 1:(ncol(df)-1))]
  fwrite(as.data.table(df), file=file)
}

#####################################################################
# Functions to create SignatureAnalyzer in-memory representations
# from standard in-memory representations.
#####################################################################


# load("SA.96.row.order.Rdata") # Defines SA.96.row.order

#' Prepare input data for SignatureAnalyzer
#'
#' @param cat96 Input catalog in in standard in-memory format
#'
#' @return Catalog in in-memory format expected by SignatureAnalyzer
#'
#' @export
SACat96 <- function(cat96) {
  stopifnot(nrow(cat96) == 96)
  SA.96.row.order <-
    c("TGTT", "TGGT", "TGCT", "TGAT", "TGTG", "TGGG", "TGCG", "TGAG",
      "TGTC", "TGGC", "TGCC", "TGAC", "TGTA", "TGGA", "TGCA", "TGAA",
      "TCTT", "TCGT", "TCCT", "TCAT", "TCTG", "TCGG", "TCCG", "TCAG",
      "TCTC", "TCGC", "TCCC", "TCAC", "TCTA", "TCGA", "TCCA", "TCAA",
      "TATT", "TAGT", "TACT", "TAAT", "TATG", "TAGG", "TACG", "TAAG",
      "TATC", "TAGC", "TACC", "TAAC", "TATA", "TAGA", "TACA", "TAAA",
      "CAAA", "CAAC", "CAAG", "CAAT", "CACA", "CACC", "CACG", "CACT",
      "CAGA", "CAGC", "CAGG", "CAGT", "CATA", "CATC", "CATG", "CATT",
      "CGAA", "CGAC", "CGAG", "CGAT", "CGCA", "CGCC", "CGCG", "CGCT",
      "CGGA", "CGGC", "CGGG", "CGGT", "CGTA", "CGTC", "CGTG", "CGTT",
      "CTAA", "CTAC", "CTAG", "CTAT", "CTCA", "CTCC", "CTCG", "CTCT",
      "CTGA", "CTGC", "CTGG", "CTGT", "CTTA", "CTTC", "CTTG", "CTTT"
    )
  rn <-rownames(cat96)
  crn <- function(rn) {
    x <-unlist(strsplit(rn, ""))
    paste(x[c(2, 4, 1, 3)], collapse = "")
  }
  nrn <- lapply(rn, crn)
  rownames(cat96) <- nrn
  cat96 <- cat96[SA.96.row.order, ]
  return(cat96)
}

#####################################################################
# Miscellaneous addtional SignatureAnalyzer-related functions
#####################################################################

#' Plot the SBS96 part of a SignatureAnalyzer COMPOSITE signature or catalog
#'
#' @param catalog Catlog or signature matrix
#'
#' @param name Name of file to print to.
#'
#' @param type See \code{\link[ICAMS]{Cat96ToPdf}}.
#'
#' @importFrom ICAMS Collapse1536To96 Cat96ToPdf
#' @export
#'
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
  Cat96ToPdf(catalog = cat96/sum(cat96), name = name, type = type)
}
