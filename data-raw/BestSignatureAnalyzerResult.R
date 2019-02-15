# Code to find the best SignatureAnalyzer result
#
#
#' Find the best SingaureAnalyzer result over a set of output directories.
#'
#' @param root.dir The directory containing the output directories.
#'
#' @result The path to the directory with the best output.
#'
#' @details As per Jaegil, we first find the most common number
#' of signatures across all the runs, and then among those
#' runs with that number of signatures, we choose the lowest
#' "evidence" (which is the negative posterior probability).
