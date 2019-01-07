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

