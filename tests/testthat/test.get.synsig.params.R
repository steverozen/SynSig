context("Getting exposure/attribution parameters and generating synthetic tumors")

test_that("synsig.params.from.attributions",
          {
            load("sa.test.param.in.Rdata")
            load("sa.test.param.out.Rdata")
            expect_equal(sa.test.param.out,
                         synsig.params.from.attributions(
                           sa.test.param.in))
          })

test_that("generate.synthetic.exposures",
          {
            load("sa.test.param.out.Rdata")
            load("sa.test.synthetic.exposures.Rdata")
            set.seed(1066)
            foo <- generate.synthetic.exposures(
              sa.test.param.out,
              num.samples = 50)
            expect_equal(foo, sa.test.synthetic.exposures)
            new.param <-
              synsig.params.from.attributions(foo)
            load("sa.test.param.in.Rdata")
            print(sa.test.param.out - new.param)
            new.delta <-
              sa.test.param.out - new.param
            load("sa.test.param.delta.Rdata")
            expect_equal(new.delta, sa.test.param.delta)
          })

