context("Test matPCs")

test_that(
  "Case matPCs.1: Missing parameters",
  {
    expect_error(matPCs())
    expect_error(matPCs(percentVar = 0.8))
  }
)
