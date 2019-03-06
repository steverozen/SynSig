context("Long test of SignatureAnalyzer4MatchedCatalogs")

test_that("SignatureAnalyzer4MatchedCatalogs", {
  # load("tmp.test.Rdata")
  load("SignatureAnalyzer4MatchedCatalogs.out.Rdata")
  set.seed(888)
  tmp.result <<-
    SignatureAnalyzer4MatchedCatalogs(
      num.runs = 2,
      signatureanalyzer.code.dir =
        "../../data-raw/SignatureAnalzyer.052418/",
      dir.root = "test.random.sigs/",
      slice = 1,
      test.only = TRUE, # Just use a small subset of the input catalog
      overwrite = TRUE,
      mc.cores = 1)
  # save(SignatureAnalyzer4MatchedCatalogs.out,
  #     file = "SignatureAnalyzer4MatchedCatalogs.out.Rdata")
  expect_equal(
    tmp.result,
    SignatureAnalyzer4MatchedCatalogs.out)
})
