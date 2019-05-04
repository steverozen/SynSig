context("Long test of SignatureAnalyzerOneRun and SignatureAnalyzer4MatchedCatalogs, > ~120 seconds on laptop")

test_that("SignatureAnalyzerOneRun", {
  # skip("This is a long test")
  load("SignatureAnalyzerOneRun.out.Rdata")
  set.seed(888)
  tmp.result <<-
    SignatureAnalyzerOneRun(
      signatureanalyzer.code.dir =
        "../../data-raw/SignatureAnalzyer.052418/",
      input.catalog = "test.random.sigs/sp.sp/ground.truth.syn.catalog.csv",
      read.catalog.function = ReadCatSNS96,
      out.dir = "test.random.sigs/sp.sp/sa.results/",
      write.signature.function = WriteCatSNS96,
      input.exposures = "test.random.sigs/sa.sa.96/ground.truth.syn.exposures.csv",
      test.only = TRUE, # Just use a small subset of the input catalog
      overwrite = TRUE)
  # save(SignatureAnalyzerOneRun.out,
  #     file = "SignatureAnalyzerOneRun.out.Rdata")

  expect_true(dir.exists("test.random.sigs/sp.sp/sa.results/"))
  expect_true(file.exists("test.random.sigs/sp.sp/sa.results/input.syn.exp.csv"))
  expect_true(file.exists("test.random.sigs/sp.sp/sa.results/sa.output.other.data.csv"))
  expect_true(file.exists("test.random.sigs/sp.sp/sa.results/sa.output.raw.exp.csv"))
  expect_true(file.exists("test.random.sigs/sp.sp/sa.results/sa.output.sigs.csv"))
  # Clean up
  dir.to.unlink <- "test.random.sigs/sp.sp/sa.results"
  unlink.ret <- unlink(dir.to.unlink, recursive = TRUE, force = TRUE)
  if (0 != unlink.ret) warning("failed to unlink", dir.to.unlink)
  expect_equal(
    tmp.result,
    SignatureAnalyzerOneRun.out)
})

test_that("SignatureAnalyzer4MatchedCatalogs", {
  # skip("This is a long test")
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

  # Clean up
  dir.to.unlink <- "test.random.sigs/sa.sa.96/sa.results"
  unlink.ret <- unlink(dir.to.unlink, recursive = TRUE, force = TRUE)
  if (0 != unlink.ret) warning("failed to unlink", dir.to.unlink)
  expect_equal(
    tmp.result,
    SignatureAnalyzer4MatchedCatalogs.out)
})
