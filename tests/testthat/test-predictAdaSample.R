context("Test predictAdaSample")

test_that(
  "Case predictAdaSample.1: Missing parameters",
  {
    cls = label = 1:12

    expect_error(predictAdaSample())
    expect_error(predictAdaSample(cls.truth = cls))
    expect_error(predictAdaSample(final = label))
  }
)
