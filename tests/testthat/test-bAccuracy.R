context("Test bAccuracy")

test_that(
  "Case bAccuracy.1: Missing parameters",
  {
    cls = label = 1:12

    expect_error(bAccuracy())
    expect_error(bAccuracy(cls.truth = cls))
    expect_error(bAccuracy(final = label))
  }
)
