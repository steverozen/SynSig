context("Mutational Signature Matching")
library(SigMisc)

test_that("Test a top level function", {
  expect_equal(
    MatchSigs2Directions(.test.extracted.sigs, .test.ground.truth.sigs),
    .test.match.out)
})
