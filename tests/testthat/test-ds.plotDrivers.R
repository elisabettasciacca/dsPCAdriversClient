# End-to-end tests for ds.plotDrivers -----------------------------------------
#
# Two layers:
#   1. Input validation (pure R, no DSLite required).
#   2. End-to-end execution via DSLite, exercising the federated PCA path
#      (dssPrincomp) and the cross-site Fisher combination.

# Input validation ------------------------------------------------------------

test_that("ds.plotDrivers rejects missing or invalid 'data' / 'vars'", {
  testthat::skip_if_not_installed("DSI")
  conn <- fake_ds_connection()

  expect_error(
    ds.plotDrivers(data = NULL, vars = "clinical", datasources = list(site1 = conn)),
    "Argument 'data' must not be NULL"
  )
  expect_error(
    ds.plotDrivers(data = "expr", vars = NULL, datasources = list(site1 = conn)),
    "Argument 'vars' must not be NULL"
  )
  expect_error(
    ds.plotDrivers(data = 1, vars = "clinical", datasources = list(site1 = conn)),
    "must be a character string"
  )
  expect_error(
    ds.plotDrivers(data = c("a", "b"), vars = "clinical",
                   datasources = list(site1 = conn)),
    "must have length 1"
  )
})

test_that("ds.plotDrivers rejects invalid datasources", {
  expect_error(
    ds.plotDrivers(data = "expr", vars = "clinical",
                   datasources = list("not a connection")),
    "must be a list of DSConnection objects"
  )
})

test_that("ds.plotDrivers rejects invalid 'type'", {
  testthat::skip_if_not_installed("DSI")
  conn <- fake_ds_connection()

  expect_error(
    ds.plotDrivers(data = "expr", vars = "clinical",
                   datasources = list(site1 = conn),
                   type = "bogus"),
    "'type' must be either 'combine' or 'split'"
  )
})


# End-to-end via DSLite -------------------------------------------------------

test_that("ds.plotDrivers (type = 'combine') returns a single ggplot across two sites", {
  testthat::skip_if_not_installed("ggplot2")

  setup <- setup_two_site_dslite()
  conns <- setup$conns
  withr::defer(DSI::datashield.logout(conns))

  p <- ds.plotDrivers(
    data        = "expr",
    vars        = "clinical",
    datasources = conns,
    type        = "combine",
    n_pc        = 3L,
    sig_cutoff  = 0.05
  )

  expect_s3_class(p, "ggplot")
  # Single (non-faceted) plot -> facet must be FacetNull
  expect_true(inherits(p$facet, "FacetNull"))
})

test_that("ds.plotDrivers (type = 'split') returns a faceted ggplot across two sites", {
  testthat::skip_if_not_installed("ggplot2")

  setup <- setup_two_site_dslite()
  conns <- setup$conns
  withr::defer(DSI::datashield.logout(conns))

  p <- ds.plotDrivers(
    data        = "expr",
    vars        = "clinical",
    datasources = conns,
    type        = "split",
    n_pc        = 3L
  )

  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$facet, "FacetWrap"))
})

test_that("ds.plotDrivers returns a results / metadata list when return_pvalues = TRUE", {

  setup <- setup_two_site_dslite()
  conns <- setup$conns
  withr::defer(DSI::datashield.logout(conns))

  out <- ds.plotDrivers(
    data           = "expr",
    vars           = "clinical",
    datasources    = conns,
    type           = "combine",
    n_pc           = 3L,
    return_pvalues = TRUE
  )

  expect_type(out, "list")
  expect_named(out, c("results", "metadata"), ignore.order = TRUE)
  expect_s3_class(out$results, "data.frame")
  expect_true(all(c("Feature", "PC", "pvalue", "Association", "Significant")
                  %in% colnames(out$results)))

  expect_equal(out$metadata$n_sites, 2L)
  expect_equal(out$metadata$meta_method, "fisher")
  # total_n must match the sum of the two simulated sites
  expect_equal(out$metadata$total_n,
               nrow(setup$site_data$site1$clinical) +
               nrow(setup$site_data$site2$clinical))

  # All p values must be in [0, 1] or NA
  expect_true(all(is.na(out$results$pvalue) |
                  (out$results$pvalue >= 0 & out$results$pvalue <= 1)))
})

test_that("ds.plotDrivers (split) returns per-site results when return_pvalues = TRUE", {

  setup <- setup_two_site_dslite()
  conns <- setup$conns
  withr::defer(DSI::datashield.logout(conns))

  out <- ds.plotDrivers(
    data           = "expr",
    vars           = "clinical",
    datasources    = conns,
    type           = "split",
    n_pc           = 3L,
    return_pvalues = TRUE
  )

  expect_type(out, "list")
  expect_named(out, c("results", "metadata"), ignore.order = TRUE)
  # One results frame per site
  expect_named(out$results, c("site1", "site2"), ignore.order = TRUE)
  expect_s3_class(out$results$site1, "data.frame")
  expect_s3_class(out$results$site2, "data.frame")

  # Site column should be present when there is more than one site
  expect_true("Site" %in% colnames(out$results$site1))
})

test_that("ds.plotDrivers (single site) forces split mode and returns a single plot", {
  testthat::skip_if_not_installed("ggplot2")

  setup <- setup_single_site_dslite()
  conns <- setup$conns
  withr::defer(DSI::datashield.logout(conns))

  # Even though type = "combine" is requested, with a single site the function
  # internally switches to "split" and returns a single (non-faceted) plot.
  p <- ds.plotDrivers(
    data        = "expr",
    vars        = "clinical",
    datasources = conns,
    type        = "combine",
    n_pc        = 3L
  )
  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$facet, "FacetNull"))

  out <- ds.plotDrivers(
    data           = "expr",
    vars           = "clinical",
    datasources    = conns,
    type           = "split",
    n_pc           = 3L,
    return_pvalues = TRUE
  )
  # Single site -> results is a single data.frame, not a named list
  expect_s3_class(out$results, "data.frame")
})

test_that("ds.plotDrivers applies p_adj after Fisher combination", {

  setup <- setup_two_site_dslite()
  conns <- setup$conns
  withr::defer(DSI::datashield.logout(conns))

  raw <- ds.plotDrivers(
    data = "expr", vars = "clinical", datasources = conns,
    type = "combine", n_pc = 3L,
    return_pvalues = TRUE
  )
  adj <- ds.plotDrivers(
    data = "expr", vars = "clinical", datasources = conns,
    type = "combine", n_pc = 3L,
    p_adj = "BH", return_pvalues = TRUE
  )

  # Same Feature x PC structure across both calls
  ok <- !is.na(raw$results$pvalue) & !is.na(adj$results$pvalue)
  expect_true(any(ok))

  # BH-adjusted p values are always >= raw p values
  expect_true(all(adj$results$pvalue[ok] >= raw$results$pvalue[ok] - 1e-12))

  # Association recomputed from adjusted p values
  expect_equal(adj$results$Association[ok],
               -log10(adj$results$pvalue[ok]),
               tolerance = 1e-10)
})

test_that("ds.plotDrivers leaves the PCA scores object on each server", {

  setup <- setup_two_site_dslite()
  conns <- setup$conns
  withr::defer(DSI::datashield.logout(conns))

  invisible(ds.plotDrivers(
    data = "expr", vars = "clinical", datasources = conns,
    type = "combine", n_pc = 3L, return_pvalues = TRUE
  ))

  symbols <- DSI::datashield.symbols(conns)
  expect_true(all(vapply(symbols, function(s) "expr_scores" %in% s,
                         logical(1))))
})
