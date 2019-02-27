context("Test summarizing extraction results in one directory")

singature.analyzer.sa.sa.COMPOSITE.out <-
  SummarizeSigOneSACOMPOSITESubdir(
    third.level.dir = "./test.sa.sa.COMPOSITE/sa.results/",
    ground.truth.exposure.name = "ground.truth.syn.exposures.csv",
    which.run = "1.run.sa.sa.COMPOSITE")

singature.analyzer.sa.sa.96.out <-
  SummarizeSigOneSA96Subdir(
    third.level.dir = "./test.sa.sa.96/sa.results/",
    ground.truth.exposure.name = "ground.truth.syn.exposures.csv",
    which.run = "1.run.sa.sa.96",
    write.png = FALSE) # ICAMS version skew / Set to false for some OS?

test_that("SummarizeSigOneSPSubdir",
          {
            load("./test.sigprofiler.sp.sp.out.Rdata")
            expect_equal(
              SummarizeSigOneSPSubdir(
                third.level.dir =
                  "./sigprofiler.sp.sp.summary.test.input/sp.results/",
                ground.truth.exposure.name = "ground.truth.syn.exposures.csv",
                write.png = FALSE), # Set this to FALSE to speed up the test.
              sigprofiler.sp.sp.out)
            if (TRUE) # Set this to FALSE to inspect the file and directory output
            {
              # Clean up
              unlink("./sigprofiler.sp.sp.summary.test.input/sp.results/summary",
                     recursive = TRUE, force = TRUE)
              unlink(
                paste0("./sigprofiler.sp.sp.summary.test.input/",
                       "sp.results/SBS96/Selected_Solution/",
                       "De_Novo_Solution/signatures.PCAWG.format.csv"))
            }
          })


