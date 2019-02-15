#' Find the best SignatureAnalyzer result over a set of output directories.
#'
#' @param root.dir The directory containing the \code{sa.run.<n>} output directories.
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
#' @return The path to the directory with the best output.
#'
#' @details As per Jaegil, we first find the most common number
#' of signatures across all the runs, and then among those
#' runs with that number of signatures, we choose the lowest
#' \code{evidence} (which is the negative posterior probability).
#'
#' @export
BestSignatureAnalyzerResult <- function(root.dir) {

  run.summaries <-
    list.files(path = root.dir,
             pattern = "sa.output.other.data.csv",
             full.names = TRUE,
             recursive = TRUE,
             include.dirs = TRUE)
  names(run.summaries) <- sub(".*(sa\\.run\\.[0-9]+).*", "\\1", run.summaries, perl = TRUE)

  process.one.summary <- function(f) {
    d <- read.csv(f, header = FALSE, row.names = 1, fill = TRUE)
    return(as.numeric(d[c("num.sigs", "evidence"), 1]))
  }

  summaries <- data.frame(lapply(run.summaries, process.one.summary))
  num.sigs <- unlist(summaries[1, ])
  evidence <- unlist(summaries[2, ])

  num.sigs.table <- data.frame(table(num.sigs))
  num.sigs.table$num.sigs <- as.numeric(as.character(num.sigs.table$num.sigs))

  # K is shorthand for num.sigs
  max.runs.for.1.K <- max(num.sigs.table$Freq)

  # There can be more than one most-common K
  most.common.Ks <- num.sigs.table[num.sigs.table$Freq == max.runs.for.1.K, "num.sigs"]

  most.common.K.indices <- which(num.sigs %in% most.common.Ks)

  num.sigs <- num.sigs[most.common.K.indices]
  evidence <- evidence[most.common.K.indices]

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
    as.numeric(sub(".*(\\d+).*", "\\1", names(num.sigs)))

  ordered.runs <- names(num.sigs)[order(names.as.numbers)]

  other.data <- run.summaries[ordered.runs[1]]
  dir.path <- sub("sa.output.other.data.csv", "", other.data, fixed = TRUE)
  return(dir.path)

}

# debug(BestSignatureAnalyzerResult)
# BestSignatureAnalyzerResult("../0syn.3.5.40.abst.v3/sa.sa.96/sa.results")

