# SigSimilarity.R
# 2019 01 06
#
# Functions to find best matches (by cosine similarity) between two
# sets of mutational signatures.
#
# Main functions of interest to users are WriteAndPlotSimilarSigs,
# which calls SigSetSimilarity.

# source("ICAMS_plotting.R")

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
#'
#' @import lsa
FindNearestSig <- function(query.sig, other.sigs) {
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
#'   is the cosine similarily between the given signature in query.sigs and
#'   the matching signature in other.sigs.
#'
#'
#' @export
#'
MatchSigs <- function(query.sigs, other.sigs) {

  stopifnot(!is.null(colnames(query.sigs)))

  Match1Sig <- function(query.sig.name) {
    query.sig <- query.sigs[ , query.sig.name]
    return(FindNearestSig(query.sig, other.sigs))
  }

  ret <-
    lapply(X = colnames(query.sigs),
           FUN = Match1Sig)
  names(ret) <- colnames(query.sigs)
  return(ret)
}

#' Calculate bidirectional closest similarities bewteen two sets of signatures
#' and the average of the similarities.
#'
#' @param sigs1 Matrix of signatures; colnames used as identifiers
#' @param sigs2 Matrix of signatures; colnames used as identifiers
#'
#' @return   #  A list with the elements:
#'    \code{avg}: the average of the cosine similarities between each signature
#'         in sigs1 and its closest match in sigs2 and the closest match
#'         between each signature in sigs2 and its closest match in sigs1
#'
#'    \code{match1}: a data frame with a signature in sigs1 in the first column,
#'            the closest match in sigs2 in the second column, and the
#'            cosine similarity between them in the third column; the rownames
#'            are the values in column 1.
#'
#'    \code{match2}: a data frame with a signature in sigs2 in the first column,
#'            the closest match in sigs1 in the second column, and the
#'            cosine similarity between them in the third column; the rownames
#'            are the values in column 1.
#'
#' \code{match1} and \code{match2} might not have the same number of rows.
#'
#' @export
#'
SigSetSimilarity <- function(sigs1, sigs2) {
  # TODO(steve): match1 and match2 can be simplified

  match1 <- MatchSigs(sigs1, sigs2)
  match2 <- MatchSigs(sigs2, sigs1)
  avg <-
    (sum(unlist(match1)) + sum(unlist(match2))) /
    (length(match1) + length(match2))

  table1 <-
    data.frame(from=names(match1),
               to=unlist(lapply(match1, names)),
               sim=unlist(match1),
               stringsAsFactors = FALSE)
  table2 <-
    data.frame(from=names(match2),
               to=unlist(lapply(match2, names)),
               sim=unlist(match2),
               stringsAsFactors = FALSE)
  return(list(avg=avg, match1=table1, match2=table2))
}

