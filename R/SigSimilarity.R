# SigSimilarity.R
#
# Functions to find best matches (by cosine similarity) between two
# sets of mutational signatures.
#
# Main functions of interest to users are WriteAndPlotSimilarSigs,
# which calls SigSetSimilarity.
#
# Will require functions from ICAMS to be passed in as arguments.

#' Find the signature in other.sigs that is nearest (by cosine similarity)
#' to query.sig.
#'
#' @param query.sig A single signature
#'
#' @param other.sigs Matrix with each column being one signature
#'
#' @return The maximum similarity between \code{query.sig} and any signature in
#' \code{other.sigs}
#'
#' @export
#' @family signature matching functions
#'
#' @import lsa
Match1Sig <- function(query.sig, other.sigs) {
  sims <-
    apply(other.sigs,
          MARGIN = 2,
          FUN = function (other.sig) {
            return(lsa::cosine(query.sig, other.sig))
            })
  max.sim <- max(sims)
  max.name <- names(sims)[sims == max.sim]
  names(max.sim) <- max.name
  return(max.sim)
}

#' Find the closest match in \code{other.sigs} for each signature in \code{query.sigs}
#'
#' @param query.sigs A signature matrix; signatures for which to find the
#'             closest match in other.sigs. The colnames are used
#'               as the identifiers of the signatures.
#'
#'
#' @param other.sigs A signature matrix; find the closest matches to
#'  a signature in this matrix.
#'  The colnames are used as the identifiers of the signatures.
#'
#'
#' @return A list with one element for each signature in query.sigs. The names
#'   of the list elements are the colnames of query.sigs. Each list element
#'   is a vector of length 1, and the name of the vector element is
#'   the name of the closest matching signature in other.sigs, and the value
#'   is the cosine similarity between the given signature in query.sigs and
#'   the matching signature in other.sigs.
#'
#' @export
#' @family signature matching functions
#'
MatchSigs1Direction <- function(query.sigs, other.sigs) {

  stopifnot(!is.null(colnames(query.sigs)))

  Match1SigInternal <- function(query.sig.name) {
    query.sig <- query.sigs[ , query.sig.name]
    return(Match1Sig(query.sig, other.sigs))
  }

  ret <-
    lapply(X = colnames(query.sigs),
           FUN = Match1SigInternal)
  names(ret) <- colnames(query.sigs)
  return(ret)
}

#' Calculate bidirectional closest similarities between two sets of signatures
#' and the average of the similarities.
#'
#' @param sigs1 Matrix of signatures; colnames are used as signature
#'   identifiers, and the colnames in \code{sigs1} should be distinguishable from
#'   those in \code{sigs2}.
#'
#' @param sigs2 Matrix of signatures; colnames are used as signature identifiers.
#'
#' @return A list with the elements:
#'
#' \code{avg}: the average of the cosine similarities between each signature in
#' \code{sigs1} and its closest match in \code{sigs2} and the closest match
#' between each signature in \code{sigs2} and its closest match in \code{sigs1}.
#'
#' \code{match1}: a data frame with rownames being signature identifiers from
#' \code{sigs1}, the signature identifier of the closest match in \code{sigs1}
#' in the 1st column, and the cosine similarity between them in the 2nd column.
#'
#' \code{match2}: a data frame with the rownames being signature identifiers
#' from \code{sigs2}, the signature identifier of the closest match in
#' \code{sigs1} in the 1st column, and the cosine similarity between them in the
#' 2nd column.
#'
#' \code{match1} and \code{match2} might not have the same number of rows.
#'
#' @export
#' @family signature matching functions
#'
MatchSigs2Directions <- function(sigs1, sigs2) {
  # TODO(Steve): match1 and match2 can be simplified

  if (is.null(colnames(sigs1))) {
    colnames(sigs1) <- paste0("Ex.", 1:ncol(sigs1))
  }

  match1 <- MatchSigs1Direction(sigs1, sigs2)
  match2 <- MatchSigs1Direction(sigs2, sigs1)
  avg <-
    (sum(unlist(match1)) + sum(unlist(match2))) /
    (length(match1) + length(match2))

  table1 <-
    data.frame(from=names(match1),
               to=unlist(lapply(match1, names)),
               sim=unlist(match1),
               stringsAsFactors = FALSE)
  rownames(table1) <- table1$from
  table1 <-table1[ , -1]

  table2 <-
    data.frame(from=names(match2),
               to=unlist(lapply(match2, names)),
               sim=unlist(match2),
               stringsAsFactors = FALSE)
  rownames(table2) <- table2$from
  table2 <-table2[ , -1]

  return(list(avg=avg, match1=table1, match2=table2))
}

#' Get the numerical parts of signature ids
#'
#' @param s A character vector
#'
#' @return A vector, each element of which is the integer
#' corresponding to the first string of digits of an element of s
NumFromId<- function(s) {
  return(
    as.numeric(
      sub("[^0123456789]*(\\d+).*", "\\1", s, perl = TRUE)))
}

#' Write a data.frame including rownames using data.table::fwrite
#'
#' @param df The data.frame to write
#' @param file Character string specifying the path to the file to create
#' @param rowname.name The colname for the first column on disk, which will
#' contain the rownames
#'
#' @export
#'
#' @import data.table
fwriteDataFrame <- function(df, file, rowname.name = "mutation.type") {
  df <- as.data.frame(df)
  df[ , rowname.name] <- rownames(df)
  df <- df[ , c(ncol(df), 1:(ncol(df)-1))]
  data.table::fwrite(as.data.table(df), file=file)
}


# There are test data generated by
#
# usethis::use_data(.test.extracted.sigs, .test.ground.truth.sigs, .test.match.out, internal = TRUE, overwrite = TRUE)
#
# usethis::use_testthat()

#' Test MatchSigs2Directions
#'
#' @details Stop on error
#' @export
#'
#' @family signature matching functions
TestMatchSigs2Directions <- function() {
  out1 <- MatchSigs2Directions(.test.extracted.sigs, .test.ground.truth.sigs)
  stopifnot(all.equal(out1, .test.match.out))
  cat("ok\n")
  invisible(out1)
}

