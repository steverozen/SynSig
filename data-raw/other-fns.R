

SigSetSimilarity <- function(sigs1, sigs2) {
  # Calculate bidirectional closest similarities bewteen two sets of signatures
  # and the average of the similarities.
  #
  # Args: sigs1 and sigs2 are matrices of sigantures; their colnames
  #       are used as idenfiers.
  #
  # Returns:
  #  A list with the elements:
  #    avg: the average of the cosine similarities between each signature
  #         in sigs1 and its closest match in sigs2 and the closest match
  #         between each signature in sigs2 and its closest match in sigs1
  #
  #    match1: a data frame with a signature in sigs1 in the first column,
  #            the closest match in sigs2 in the second column, and the
  #            cosine similarity between them in the third column; the rownames
  #            are the values in column 1. TODO(steve): can be simplified
  #            by dropping column 1.
  #
  #    match2: a data frame with a signature in sigs2 in the first column,
  #            the closest match in sigs1 in the second column, and the
  #            cosine similarity between them in the third column; the rownames
  #            are the values in column 1.
  #
  #    In general, match1 and match2 will not have the same number of rows.

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

# Get the numerical part of signature ids, e.g. c("SBS3", "SBS10") -> c(3, 10)
# (Actually, just gets the numbers corresponding to the first stretch of digits
# in each element of x)
NumFromId<- function(s) {
  return(
    as.numeric(
      sub("[^0123456789]*(\\d+).*", "\\1", s, perl = TRUE)))
}

# TODO(steve): move this to functionality to ICAMS plotting
# Call Cat96ToPdf, with the argument "id" set to the column names of the catalog
Cat96ToPdf0 <- function(catalog, name, type = "signature") {
  invisible(
    Cat96ToPdf(
      catalog = catalog,
      name = name,
      type = type, id = colnames(catalog)))
}

WriteAndPlotSimilarSigs <-
  function(ex.sigs, gt.sigs, plot.fn, file.prefix = '') {
    # Args:
    #   ex.sigs: Newly extracted signatures to be compared to gt.sigs
    #            (actually, this is more general)
    #
    #   gt.sigs: "ground truth" signatures
    #
    #   plot.fn: Function to plot signatures. For "typical" 96-channel
    #            signatures use Cat96ToPdf0
    #
    # Returns a list with the elments:
    #   avg, match1, match2: as for SigSetSimilarity, with sigs1 being the
    #                        extracted signatures and sigs2 being the
    #                        ground truth signatures (gt.sig)
    #
    #   ex.sigs: A copy of the input
    #
    #   gt.sigs: A copy of the input
    #
    #   file.prefix: A prefix string for output files, can be used to put all
    #                output files in one directory, e.g. by setting file.prefix
    #                to "my_directory/something".
    #
    #
    # Side effects:
    #
    #   1. Plots the extracted signatures with identifiers showing the nearest
    #   ground-truth signatures in a location defined by file.prefix
    #
    #   2. Saves the match tables to .csv files and overview to .csv or other
    #   text file.

    if (is.null(colnames(ex.sigs))) {
      colnames(ex.sigs) <- paste0("Ex.", 1:ncol(ex.sigs))
    }
    sim <- SigSetSimilarity(ex.sigs, gt.sigs)

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
    # nums <-
    #  as.numeric(sub("[^\\d]*(\\d+).*", "\\1",
    #                 sim$match1$to, perl = T))
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
