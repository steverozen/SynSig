#' diff new directory / files against regression data for testing.
#'
#' @param dirname the root name of the directories to diff.
#'
#' @return The return value of diff command.

Diff4SynDataSets <- function(dirname) {
  regressdirname <- paste0("long.test.regression.data/", dirname)
  if (!file.exists(regressdirname)) stop(regressdirname, " does not exist")
  tmpdirname <- paste0("tmp.", dirname)
  if (!file.exists(tmpdirname)) stop(tmpdirname, " does not exist")
  cmd <- paste("diff -rq", tmpdirname, regressdirname)
  cmd.result <- system(cmd)
  if (cmd.result == 0) {
    # No differences
    unlink.res <- unlink(tmpdirname, recursive = TRUE, force = TRUE)
    if (unlink.res != 0) warning("failed to unlink ", tmpdirname)
  }
  return(cmd.result)
}
