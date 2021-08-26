context("Test multiAdaSampling")

test_that(
  "Case multiAdaSampling.1: Missing parameters",
  {
    dat = matrix(rnorm(1200), nrow= 12)
    label = 1:ncol(dat)
    expect_error(multiAdaSampling())
    expect_error(multiAdaSampling(data = dat))
    expect_error(multiAdaSampling(label = label))
  }
)
