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

#' Write a SignatureAnalyzer catalog or signature matrix to disk
#'
#' @param ct A catalog or signature matrix
#'
#' @param path Path to file to write
#'
#' @importFrom utils write.csv
#'
#' @export
#'
WriteCatCOMPOSITE <- function(ct, path) {
  write.csv(ct, file = path, row.names = TRUE,
            quote = FALSE)
}

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
#' @param catalog Catalog or signature matrix
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

#' Standardize SignatureAnalyzer signature names
#'
#' For example, change \code{BI_COMPOSITE_SNV_SBS83_P}
#' to \code{BI_COMPOSITE_SBS83_P}
#'
#' @param sig.names Vector of signature names
#'
#' @return Vector of signatures names with "_SNV" removed.
#'
#' @export
FixSASigNames <- function(sig.names) {
  return(gsub("_SNV_", "_", sig.names, fixed = TRUE))
}

# Add this doc to the signature name fixup for SA:
#
# . This is necessary because
#' for COMPOSITE signatures we rbind coordinated
#' "SNV", "DNP", and "INDEL" signatures. See
#' \code{\link{ReadSASigsCOMPOSITE}}.


#' Function for running SignatureAnalyzer on a file containing
#' a catalog.
#'
#' It is necessary to setwd() to a folder containing
#' SignatureAnalyzer (right now SignatureAnal\strong{zy}er.052418,
#' note transposed "zy")
#' and then call
#' \preformatted{
#' options( warn = 0 )
#' here <- getwd()
#' setwd("SignatureAnalzyer.052418/") # Or the appropriate directory
#' INPUT <- "INPUT_SignatureAnalyzer/"
#' source("SignatureAnalyzer.PCAWG.function.R")
#' setwd(here) # This is necessary because the caller
#'             # as specified input and output locations
#'             # realtive to here.
#' RunSignatureAnalyzerOnSyn(...)
#' }
#'
#' TODO(Steve): see if we can do the source of SA inside
#' this function.
#'
#' @param input.catalog File containing input catalog.  Columns are
#' samples (tumors), rows are signatures.  SignatureAnalyzer does
#' not care about the row names (I think) TODO(Steve): check this.
#'
#' @param read.catalog.function Function taking a file path as
#' its only argument and returning a catalog as a numeric matrix.
#'
#' @param out.dir Directory that will be created for the output;
#' abort if it already exits.  Log files will be in
#' \code{paste0(out.dir, "/tmp")}.
#'
#' @param write.signature.function Function with first argument the
#' signatures generated by SignatureAnalyzer and second argument
#' the file to write to.
#'
#' @param maxK The maximum number of signatures to consider
#' extracting.
#'
#' @param tol Controls when SignatureAnalyzer will terminate
#' its search; \code{tol} was 1.e-05 for the PCAWG7 analysis.
#'
#' @param test.only If true, only analyze the first 10 columns
#' read in from \code{input.catalog}.
#'
#' @return The output of SignatureAnalyzer, invisibly.
#'
#' @details Creates the file \code{paste0(out.dir, "/sa.output.rdata")}.
#'
#' TODO(Steve): review the output format from SignatureAnalyzer -- ? names of sigs just come from the input?
#'
#' @export
#'
#' @importFrom utils write.csv
RunSignatureAnalyzerOnFile <-
  function(input.catalog,
           read.catalog.function,
           out.dir,
           write.signature.function,
           maxK = 30, tol = 1e-7, test.only = FALSE) {

    syn.data <- read.catalog.function(input.catalog)

    if (test.only) syn.data <- syn.data[ , 1:10]

    if (dir.exists(out.dir)) stop(out.dir, "already exits")
    dir.create(out.dir)
    # TEMPORARY is a global required by SignatureAnalyzer
    TEMPORARY <<- paste0(out.dir, "/tmp/")
    dir.create(TEMPORARY)

    # BayesNMF.L1W.L2H is defined by the statement
    # source("SignatureAnalyzer.PCAWG.function.R") above.
    out.data <-
      BayesNMF.L1W.L2H(syn.data, 200000, 10, 5, tol, maxK, maxK, 1)

    sigs <- out.data[[1]]
    sigs <- sigs[   , colSums(sigs) > 0]

    exp <- out.data[[2]]
    exp <- exp[rowSums(exp) > 0, ]

    save(out.data, sigs,
         file = paste0(out.dir, "/sa.output.rdata"))

    write.signature.function(sigs,
                             paste0(out.dir, "/sa.output.sigs.csv"),
                             strict = TRUE)

    WriteExposure(exp, file = paste0(out.dir, "/sa.output.exp.csv"))
    # write.csv(exp, file = paste0(out.dir, "/sa.output.exp.csv"))

    invisible(out.data)
  }
# out.data[[1]] is the signatures
# out.data[[2]] is the exposures
# [[3]] is numeric, length 1, e.g. 56.6531
# [[4]] -5884.391
# 5, 30 numbers (length same as maxK?)
# 6, 1 number, 10545.61

