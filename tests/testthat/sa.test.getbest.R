context("Getting best result from replicate SignatureAnalyzer runs")

test_that(
  "Test 1 for BestSignatureAnalyzerResult",
  expect_equivalent(
    BestSignatureAnalyzerResult("test.sa.results1/", verbose = TRUE),
    "test.sa.results1//sa.run.12/" )
)

test_that(
  "Test 2 for BestSignatureAnalyzerResult",
  expect_equivalent(
    BestSignatureAnalyzerResult("test.sa.results2/", verbose = TRUE),
    "test.sa.results2//sa.run.12/"
  )
)

test_that(
  "Test 3 for BestSignatureAnalyzerResult",
  expect_equivalent(
    BestSignatureAnalyzerResult("test.sa.results3/", verbose = TRUE),
    "test.sa.results3//sa.run.3/"
  )
)


#' Function 4 to test BestSignatureAnalyzerResult
#' Should return path to sa.run.1
Besttest4 <- function() {
  debug(BestSignatureAnalyzerResult)
  best.dir <-
    BestSignatureAnalyzerResult("data-raw/test.sa.results4/", verbose = TRUE)
  cat("\n", best.dir, "\n")
}
