# SigSimilarity.R
# 2019 01 07
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
#'   is the cosine similarily between the given signature in query.sigs and
#'   the matching signature in other.sigs.
#'
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
#' @family signature matching functions
#'
MatchSigs2Directions <- function(sigs1, sigs2) {
  # TODO(steve): match1 and match2 can be simplified

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
  table2 <-
    data.frame(from=names(match2),
               to=unlist(lapply(match2, names)),
               sim=unlist(match2),
               stringsAsFactors = FALSE)
  return(list(avg=avg, match1=table1, match2=table2))
}


# , e.g. c("SBS3", "SBS10") -> c(3, 10)
# (Actually, just gets the numbers corresponding to the first stretch of digits
# in each element of x)
#' Get the numerical parts of signature ids
#'
#' @param s A character vector
#'
#' @return A vector, each element of which is the integer
#' corresponding to the first string of digits of an element of s
#'
#' @example NumFromId(c("AB12c23", "AB10"))
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


#' Run \code{MatchSigs2Directions}, then
#' plot its results and write them as .csv files.
#'
#' @param ex.sigs Newly extracted signatures to be compared to gt.sigs
#            (actually, this is more general)
#
#' @param gt.sigs "Ground truth" signatures
#'
#' @param plot.fn Function to plot signatures. For "typical" 96-channel
#            signatures use \code{Cat96ToPdf0}
#
#' @param file.prefix A prefix string for output files, can be used to put all
#                output files in one directory, e.g. by setting file.prefix
#                to "my_directory/something"
#'
#' @return A list with the elments \code{avg}, \code{match1},
#' \code{match2} as for \code{SigSetSimilarity}, with \code{match1}
#'  being matches for the the extracted signatures (\code{ex.sigs})
#'  and \code{match2} being the
#'  matches for the ground truth signatures (\code{gt.sigs}). The return list
#'  also echos the input arguments \code{ex.sigs} and \code{gt.sigs}.
#
#' @export
#' @family signature matching functions
#'
#' @details Uses \code{plot.fn} to
#' plot the extracted signatures with identifiers showing the
#' nearest ground-truth signatures in a location defined by
#' \code{file.prefix}. Saves the match tables to
#' .csv file with suffixes \code{match.extracted.to.gt.csv}
#' and \code{match.gt.to.extracted.csv}.
#'
MatchSigsThenWriteAndPlot <-
  function(ex.sigs, gt.sigs, plot.fn, file.prefix = '') {

    if (is.null(colnames(ex.sigs))) {
      colnames(ex.sigs) <- paste0("Ex.", 1:ncol(ex.sigs))
    }
    sim <- MatchSigs2Directions(ex.sigs, gt.sigs)

    # Write the "match" tables as .csv files
    fwriteDataFrame(sim$match1,
                    paste0(file.prefix, "match.extracted.to.gt.csv"),
                    rowname.name = "Extracted.sig")
    fwriteDataFrame(sim$match2,
                    paste0(file.prefix, "match.gt.to.extracted.csv"),
                    rowname.name = "Ground.truth.sig")

    # TODO(steve) Document the complexity below; mostly it deals
    # with setting up plotting that is easy(?) to interpret.
    labels <- character(ncol(ex.sigs))
    names(labels) <- colnames(ex.sigs)
    nums <- NumFromId(sim$match1$to)
    reordered.ex <- colnames(ex.sigs)[order(nums)]
    ex.sigs.x <- ex.sigs[ , order(nums)]
    bestmatch.id <- sim$match1[reordered.ex, "to"]
    bestmatch.sim <- sim$match1[reordered.ex, "sim"]
    bestmatch.sim <- round(bestmatch.sim, digits=4)
    init.labels <-
      paste0(reordered.ex, " (", bestmatch.id, " ", bestmatch.sim, ")")
    names(init.labels) <- reordered.ex
    laggards <- setdiff(sim$match2$from, bestmatch.id)
    # Falling back to a loop here:
    for (lag in laggards) {
      my.ex.id  <- sim$match2[lag, "to"]
      my.ex.sim <- round(sim$match2[lag, "sim"], digits = 4)
      init.labels[my.ex.id] <-
        paste0(init.labels[my.ex.id],
               " (", lag, " ", my.ex.sim, ")")
    }
    colnames(ex.sigs.x) <- init.labels
    plot.file <- paste0(file.prefix, "extracted.sigs.pdf")
    # cat("Plotting to", plot.file, "\n")
    plot.fn(ex.sigs.x, plot.file)

    sim$ex.sigs <- ex.sigs
    sim$gt.sigs <- gt.sigs
    invisible(sim)
  }

# There are test data generated by
# usethis::use_data(.test.extracted.sigs, .test.ground.truth.sigs,
#                   internal = TRUE, overwrite = TRUE)
