context("Generating synthetic tumors")

test_that("synsig.params.from.attributions",
          {
            cat("HERE I AM", getwd(), "\n")
            load("sa.test.param.in.Rdata")
            load("sa.test.param.out.Rdata")
            expect_equal(sa.test.param.out,
                         synsig.params.from.attributions(sa.test.param.in))
          })
