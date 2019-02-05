
#' @param plot.fn Function to plot signatures. For "typical" 96-channel
#            signatures use \code{Cat96ToPdf0}
#
#' @param file.prefix A prefix string for output files, can be used to put all
#                output files in one directory, e.g. by setting file.prefix
#                to "my_directory/something"
#'
#'
#'
#'@details Uses \code{plot.fn} to
#' plot the extracted signatures with identifiers showing the
#' nearest ground-truth signatures in a location defined by
#' \code{file.prefix}. Saves the match tables to a
#' .csv file with suffixes \code{match.extracted.to.gt.csv}
#' and \code{match.gt.to.extracted.csv}. \strong{IMPORTANT:}
#' This function does not check that all the \code{gt.sigs}
#' were "exposed" in generating the synthetic data from
#' which \code{ex.sigs} were generated. TODO(Steve): Possibly
#' add this check.
#'
#'#'
#'
#' @importFrom utils write.csv
#'
#'

plot.file <- paste0(file.prefix, "extracted.sigs.pdf")
# cat("Plotting to", plot.file, "\n")
plot.fn(ex.sigs.x, plot.file, type = "signature")
# Write the "match" tables as .csv files
write.csv(sim$match1,
          paste0(file.prefix, "match.extracted.to.gt.csv"),
          row.names = TRUE)
write.csv(sim$match2,
          paste0(file.prefix, "match.gt.to.extracted.csv"),
          row.names = TRUE)

