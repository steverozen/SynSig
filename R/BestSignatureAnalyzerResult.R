#' Find the best SignatureAnalyzer result over a set of output directories.
#'
#' @param sa.results.dir The directory containing the \code{sa.run.<n>} output directories.
#'  For example, contents of \code{../syn.3.5.40.abst.v3/sa.sa.96/sa.results/}
#'  must be the directories
#'  \preformatted{
#'  sa.run.1
#'  sa.run.2
#'  sa.run.3
#'  ...
#'  }
#' each containing
#' \preformatted{
#' sa.output.other.data.csv
#' sa.output.sigs.csv
#' ...
#' }
#'
#' (This function only looks for subdirectories with names of the form
#' \code{sa.run.<n>}.)
#'
#' Then
#'  \preformatted{
#'  #'  BestSignatureAnalyzerResult(
#'    "../syn.3.5.40.abst.v3/sa.sa.96/sa.results")
#'  }
#'  returns the path to the sa.run.<n> directory with
#'  the output from the best SignatureAnalyzer run, e.g.
#' \preformatted{"../syn.3.5.40.abst.v3/sa.sa.96/sa.results/sa.run.3/"}
#'
#' @param verbose If TRUE print informative messages to standard output.
#'
#' @return The path to the directory with the best output as a
#' string, with the list directories examined as the attribute
#' \code{run.directories}.
#'
#' @details As per Jaegil, we first find the most common number
#' of signatures across all the runs, and then among those
#' runs with that number of signatures, we choose the lowest
#' \code{evidence} (which is the negative posterior probability).
#'
#' @keywords internal
BestSignatureAnalyzerResult <- function(sa.results.dir,
                                        verbose = FALSE) {
  me <- match.call()[[1]]

  run.directories <-
    list.files(path = sa.results.dir,
               pattern = "run",
               full.names = TRUE,
               recursive = FALSE,
               include.dirs = TRUE)
  if (length(run.directories) == 0) {
    stop("No *run* directories in", sa.results.dir)
  }

  run.summaries <-
    list.files(path = sa.results.dir,
             pattern = "sa.output.other.data.csv",
             full.names = TRUE,
             recursive = TRUE,
             include.dirs = TRUE)
  run.summaries <-
    grep("sa.run.\\d+", run.summaries, value = TRUE, perl = TRUE)
  names(run.summaries) <-
    sub(".*(sa\\.run\\.[0-9]+).*", "\\1", run.summaries, perl = TRUE)

  process.one.summary <- function(f) {
    d <- read.csv(f, header = FALSE, row.names = 1, fill = TRUE)
    return(as.numeric(d[c("num.sigs", "evidence"), 1]))
  }

  summaries <- data.frame(lapply(run.summaries, process.one.summary))
  num.sigs <- unlist(summaries[1, ])
  if (verbose) {
    cat("\n", me, "vector of K:\n")
    print(num.sigs)
    cat("\n\n")}
  evidence <- unlist(summaries[2, ])
  if (verbose) {
    cat("\n", me, "vector of evidence:\n")
    print(evidence)
    cat("\n\n")}

  num.sigs.table <- data.frame(table(num.sigs))
  num.sigs.table$num.sigs <- as.numeric(as.character(num.sigs.table$num.sigs))

  # K is shorthand for num.sigs
  max.runs.for.1.K <- max(num.sigs.table$Freq)

  # There can be more than one most-common K
  most.common.Ks <- num.sigs.table[num.sigs.table$Freq == max.runs.for.1.K, "num.sigs"]
  if (verbose)
    cat("\n", me,
        " most common numbers of extracted signatures:",
        most.common.Ks, "\n\n")

  most.common.K.indices <- which(num.sigs %in% most.common.Ks)

  num.sigs <- num.sigs[most.common.K.indices]
  evidence <- evidence[most.common.K.indices]
  if (verbose) {
    cat("\n", me, "evidences for most common Ks:\n")
    print(evidence)
    cat("\n\n")
  }

  best.evidence.indices <- which(evidence == min(evidence))
  num.sigs <- num.sigs[best.evidence.indices]
  evidence <- evidence[best.evidence.indices]

  # There still might be mutliple runs with equally good evidence,
  # so we take the ones with the most signatures.

  max.K <- max(num.sigs)
  max.K.indices <- which(num.sigs == max.K)

  num.sigs <- num.sigs[max.K.indices]
  evidence <- evidence[max.K.indices]

  # There still might be equally good runs, so we take the
  # lowest numbered run.

  names.as.numbers <-
    as.numeric(sub("[^\\d]*(\\d+).*", "\\1", names(num.sigs),
                   perl = TRUE))

  ordered.runs <- names(num.sigs)[order(names.as.numbers)]

  other.data <- run.summaries[ordered.runs[1]]
  dir.path <- sub("sa.output.other.data.csv", "", other.data, fixed = TRUE)

  # Do not want to change the interface, but this information
  # is useful downstream
  attr(dir.path, "run.directories") <- run.directories

  return(dir.path)

}

#' Find best result and make a copy of the folder
#'
#' @param sa.results.dir See \code{\link{BestSignatureAnalyzerResult}}
#'
#' @param verbose See \code{\link{BestSignatureAnalyzerResult}}
#'
#' @param overwrite If TRUE overwrite existing "best.run"
#'
#' @return The path of the best directory that was copied as a
#' string, with the list directories examined as the attribute
#' \code{run.directories}.
#'
#' @export
#'

CopyBestSignatureAnalyzerResult <-
  function(sa.results.dir,
           verbose = FALSE,
           overwrite = FALSE) {
    best <- BestSignatureAnalyzerResult(sa.results.dir, verbose)
    target.dir <- paste0(sa.results.dir, "/best.run")
     # file.copy.ret <- file.copy(
    #  from = best, to = target.dir, overwrite = TRUE)
    #copyDirectory(
    #  from = best, to = target.dir, overwrite = overwrite, recursive = TRUE)
    if (dir.exists(target.dir)) {
      if (!overwrite) stop("Directory", target.dir, " exists and !overwrite")
    } else {
      create.result <- dir.create(target.dir)
      if (!create.result) warning("problem creating", target.dir)
    }
    best.exp <- paste0(best, "sa.output.exp.csv")
    if (!file.exists(best.exp)) stop(best.exp, " does not exist")
    file.copy(from = best.exp, to = target.dir, overwrite = overwrite)
    file.copy(from = paste0(best, "sa.output.sigs.csv"),
              to = target.dir, overwrite = overwrite)
    file.copy(from = paste0(best, "sa.output.other.data.csv"),
              to = target.dir, overwrite = overwrite)

    info.file.name <- paste0(target.dir, "/info.txt")
    cat("This is a copy of", best, "\n\n", file = info.file.name)
    run.directories <- attr(best, "run.directories")
    cat("Based on", length(run.directories), "runs:\n",
        file = info.file.name, append = TRUE)
    cat(run.directories, "\n", file = info.file.name, append = TRUE)

    return(best)
  }
