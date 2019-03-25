#' diff new directory / files against regression data for testing.
#'
#' @param dirname the root name of the directories to diff.
#'
#' @param unlink if \code{TRUE} unlink \code{tmpdirname}, but do not unlink
#' if there are diffs.
#'
#' @return Invisibly, the return value of the diff command.
#'
#' @export

Diff4SynDataSets <- function(dirname, unlink) {
  regressdirname <- paste0("long.test.regression.data/", dirname)
  if (!file.exists(regressdirname)) stop(regressdirname, " does not exist")
  tmpdirname <- paste0("tmp.", dirname)
  if (!file.exists(tmpdirname)) stop(tmpdirname, " does not exist")
  cmd.result <- system2("diff", c("-rq", tmpdirname, regressdirname))
  if (cmd.result == 0) {
    # No differences
    if (unlink) {
      unlink.res <- unlink(tmpdirname, recursive = TRUE, force = TRUE)
      if (unlink.res != 0) warning("failed to unlink ", tmpdirname)
    }
  }
  invisible(cmd.result)
}
