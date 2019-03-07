context("Test generation of ranodom signatures")

test_that("CreateRandomMutSigProfiles", {
  load("rand.syn.96.sigs.Rdata")
  set.seed(5)
  expect_equal(
    CreateRandomMutSigProfiles(
      ICAMS::catalog.row.order[["SNS96"]], 5, "prefix"),
    rand.syn.96.sigs)
})
