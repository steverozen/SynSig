#####################################################################
# Read SignatureAnalyzer signature catalogs from the
# PCAWG7 format (as found on Synapse) into standard
# in-memory representations.
#
# Used to generate the package variables and not tested since.
#####################################################################

library(ICAMS)
library(SynSig) # ReadSASigCOMPOSITE, ReadSASig96
library(usethis)

library(data.table)

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

#' Read signatures in SignatureAnalyzer format as found on
#' the PCAWG7 Synapse web site.
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

ReadAnySASig <- function(path, canonicalize.rows, exp.nrow, row.order) {
  dt <- as.data.frame(fread(path))
  stopifnot(nrow(dt) == exp.nrow)
  feature <- which(names(dt) == "feature")
  feature.names <- dt[ , feature]
  dt <- dt[ , -1 * feature]

  rownames(dt) <- canonicalize.rows(feature.names)

  out <- dt[row.order, ]

  return(as.matrix(out))
}

#' Read SignatureAnalyzer signatures from PCAWG7 / Synapse.
#'
#' There is no custom catalog reading function for SignatureAnalyzer.
#'
#' \code{ReadSASig96} Read a 96 SNS catalog from path
#'
#' \code{ReadSASig1536} Read a 1536 SNS catalog from path
#'
#' \code{ReadSASigDBS} Read a 78 DNS catalog from path
#'
#' \code{ReadSASigID} Read a ID (insertion/deletion) catalog from path
#' Please take note that deletion Repeat Size ranges from 0 to 5+
#' in the catalog, but for plotting and end user documentation
#' it ranges from 1 to 6+.
#'
#' @param path The file to read from.
#'
#' @seealso \code{\link{ReadSASigCOMPOSITE}}
#'
#' @return A catalog in canonical in-memory format.
#'
#' @name ReadSASig
NULL

#' @rdname ReadSASig
#' @export
ReadSASig96 <- function(path) {
  ReadAnySASig(path, CanonicalizeSAx, 96, ICAMS::catalog.row.order[["SNS96"]])
}

#' @rdname ReadSASig
#' @export
ReadSASig1536 <- function(path) {
  ReadAnySASig(path, CanonicalizeSAx, 1536, ICAMS::catalog.row.order[["SNS1536"]])
}

#' @rdname ReadSASig
#' @export
ReadSASigDBS <- function(path) {
  ReadAnySASig(path, CanonicalizeSADBS, 78, ICAMS::catalog.row.order[["DNS78"]])
}

#' @rdname ReadSASig
#' @export
ReadSASigID <- function(path) {
  ReadAnySASig(path, CanonicalizeSAID, 83, ICAMS::catalog.row.order[["ID"]])
}

#' Read SignatureAnalyzer COMPOSITE signatures from fixed files
#'
#' This function was used in loading SignatureAnalyzer
#' signatures from PCAWG7 / Synapse, and has not been
#' tested since.
#'
#' @return A matrix (data.frame?) of COMPOSITE signatures
#'
#' @export
ReadSASigCOMPOSITE <- function() {

  sa1536.file <-
    "SignatureAnalyzer_COMPOSITE_SBS_W1536.signature.031918.txt"
  saDBS.file  <-
    "SignatureAnalyzer_COMPOSITE_DBS.signature.042018.txt"
  saID.file   <-
    "SignatureAnalyzer_COMPOSITE_ID.signature.031918.txt"

  sa1536.sig  <- ReadSASig1536(sa1536.file)
  saDBS.sig   <- ReadSASigDBS(saDBS.file)
  saID.sig    <- ReadSASigID(saID.file)

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


# COMPOSITE signature profiles
sa.COMPOSITE.sigs <- ReadSASigCOMPOSITE()


# 96-channel signature profiles
sa.96.sigs <-
  ReadSASig96("SignatureAnalyzer_COMPOSITE_SBS_W96.signature.031918.txt")
colnames(sa.96.sigs) <- FixSASigNames(colnames(sa.96.sigs))

usethis::use_data(sa.COMPOSITE.sigs, sa.96.sigs)
