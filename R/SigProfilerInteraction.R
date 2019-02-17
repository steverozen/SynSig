# SigProfilerInteraction.R

#' Read a file containing signatures extracted by SigProfiler/Python
#'
#' @param file The name of the file to read.
#'
#' @return The corresponding signature matrix in standard internal
#' representation.
#'
#' @importFrom utils read.table
#'
#' @export

ReadSigProfilerSig96 <- function(file) {
  x <- read.table(file, sep = "\t", as.is = TRUE, header = TRUE)
  n <- x[ ,1]
  x <- x[ , -1]

  # A[C>T]A --> ACAT
  new.n <-
    paste0(substr(n, 1, 1), substr(n, 3, 3), substr(n, 7, 7), substr(n, 5, 5))

  rownames(x) <- new.n

  x <- x[ICAMS:::.catalog.row.order96, ]

  return(x)

}

## Turn this into a test
## ReadSigProfilerSig96("example-SP-signatures.txt")



