context("Test create synthetic catalogs representing several tumor types")

test_that("CreateMixedTumorTypeSyntheticData", {
  load("mixed.types.Rdata")
  set.seed(191906)
  num.syn.tumors <- 10
  cancer.types <- c("Bladder-TCC", "Skin-Melanoma")
  expect_equal(
    CreateMixedTumorTypeSyntheticData(
      top.level.dir = "./Bladder-Melanoma-test",
      cancer.type.strings = cancer.types,
      num.syn.tumors = num.syn.tumors,
      overwrite = TRUE
    ),
    mixed.types)

    if (TRUE) # Set this to FALSE to inspect the file and directory output
  {
    # Clean up
    unlink("./Bladder-Melanoma-test",
           recursive = TRUE, force = TRUE)
  }
})
