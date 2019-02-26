context("Getting exposure/attribution parameters and generating synthetic tumors")

test_that("GetSynSigParamsFromExposures",
          {
            load("sa.test.param.in.Rdata")
            load("sa.test.param.out.Rdata")
            expect_equal(sa.test.param.out,
                         GetSynSigParamsFromExposures(
                           sa.test.param.in))
          })

test_that("GenerateSyntheticExposures",
          {
            load("sa.test.param.out.Rdata")
            load("sa.test.synthetic.exposures.Rdata")
            set.seed(1066)
            foo <- GenerateSyntheticExposures(
              sa.test.param.out,
              num.samples = 50)
            expect_equal(foo, sa.test.synthetic.exposures)
            new.param <-
              GetSynSigParamsFromExposures(foo)
            load("sa.test.param.in.Rdata")
            new.delta <-
              sa.test.param.out - new.param
            # print(new.delta)
            load("sa.test.param.delta.Rdata")
            expect_equal(new.delta, sa.test.param.delta)
          })

