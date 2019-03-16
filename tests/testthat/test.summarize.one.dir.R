context("Test summarizing extraction results in one directory")

test_that("signature.analyzer.sa.sa.COMPOSITE.out", {
  load("signature.analyzer.sa.sa.COMPOSITE.out")
  expect_equal(
    SummarizeSigOneSACOMPOSITESubdir(
      third.level.dir = "./test.sa.sa.COMPOSITE/sa.results/",
      ground.truth.exposure.name = "ground.truth.syn.exposures.csv",
      which.run = "1.run.sa.sa.COMPOSITE",
      overwrite = FALSE),
    signature.analyzer.sa.sa.COMPOSITE.out)
  if (TRUE) # Set this to FALSE to inspect the file and directory output
  {
    # Clean up
    unlink("./test.sa.sa.COMPOSITE/sa.results/summary",
           recursive = TRUE, force = TRUE)
  }
})

test_that("SummarizeSigOneSA96Subdir", {
  load("signature.analyzer.sa.sa.96.out.Rdata")
  expect_equal(

    SummarizeSigOneSA96Subdir(
      third.level.dir = "./test.sa.sa.96/sa.results/",
      ground.truth.exposure.name = "ground.truth.syn.exposures.csv",
      which.run = "1.run.sa.sa.96",
      write.png = FALSE,  # Set to FALSE for some operating systems
      overwrite = TRUE),
    signature.analyzer.sa.sa.96.out)
  if (TRUE) # Set this to FALSE to inspect the file and directory output
  {
    # Clean up
    file.to.unlink <- "./test.sa.sa.96/sa.results/summary"
    if (!file.exists(file.to.unlink)) cat("Coding error,", file.to.unlink,
                                          "does not exist\n")
    res <- unlink(file.to.unlink,
           recursive = TRUE, force = TRUE)
    if (res != 0) cat("Failed to unlink", file.to.unlink, "\n")
  }
})

test_that("SummarizeSigOneSPSubdir", {
  load("./test.sigprofiler.sp.sp.out.Rdata")
  expect_equal(
    SummarizeSigOneSPSubdir(
      third.level.dir =
        "./sp.sum.test.in/sp.results/", overwrite = T,
      ground.truth.exposure.name = "ground.truth.syn.exposures.csv",
      write.png = FALSE), # Set this to FALSE to speed up the test.
    sigprofiler.sp.sp.out)
  if (TRUE) # Set this to FALSE to inspect the file and directory output
  {
    # Clean up
    file1 <- "./sp.sum.test.in/sp.results/summary"
    if (!file.exists(file1)) cat("coding error, wrong file\n")
    res <- unlink(file1, recursive = TRUE, force = TRUE)
    if (res != 0) cat("Failed to unlink", file1, "\n")
    res <- unlink(
      paste0("./sp.sum.test.in/",
             "sp.results/SBS96/Selected_Solution/",
             "De_Novo_Solution/signatures.PCAWG.format.csv"))
    if (res != 0) cat("Failed to unlink several files\n")
  }
})


