# Note: Full testing requires DataSHIELD server connection
# For unit tests, mock DSI functions

test_that("ds.plotDrivers validates inputs", {
  expect_error(
    ds.plotDrivers(pcs = NULL, vars = "clinical"),
    "pcs must be a single character string"
  )

  expect_error(
    ds.plotDrivers(pcs = "pca", vars = NULL),
    "vars must be a single character string"
  )
})
