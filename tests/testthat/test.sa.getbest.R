context("Getting best result from replicate SignatureAnalyzer runs")

test_that(
  "Test 1 for BestSignatureAnalyzerResult",
  expect_equivalent(
    BestSignatureAnalyzerResult("test.sa.results1/"),
    "test.sa.results1//sa.run.12/" )
)

test_that(
  "Test 2 for CopyBestSignatureAnalyzerResult",
  {
    expect_equivalent(
      CopyBestSignatureAnalyzerResult("test.sa.results1/"),
      "test.sa.results1//sa.run.12/" )
    # Clean up. Otherwise the copy will fail on the
    # next run of the test.
    unlink("test.sa.results1/best.run", recursive = TRUE)
  }
)

test_that(
  "Test 3 for BestSignatureAnalyzerResult",
  expect_equivalent(
    BestSignatureAnalyzerResult("test.sa.results2/"),
    "test.sa.results2//sa.run.12/"
  )
)

test_that(
  "Test 4 for BestSignatureAnalyzerResult",
  expect_equivalent(
    BestSignatureAnalyzerResult("test.sa.results3/"),
    "test.sa.results3//sa.run.3/"
  )
)

test_that(
  "Test 5 for BestSignatureAnalyzerResult",
  expect_equivalent(
    BestSignatureAnalyzerResult("test.sa.results4/"),
    "test.sa.results4//sa.run.1/"
  )
)
