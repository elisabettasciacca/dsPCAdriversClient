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
  expect_true(all(c("Feature", "PC", "pvalue", "Significant")
                  %in% colnames(out$results)))
  # No p_adj was requested, so there is no pvalue_adj column.
  expect_false("pvalue_adj" %in% colnames(out$results))

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

  # pvalue_adj is added only when p_adj is requested - not a spurious extra
  # or missing column either way.
  expect_named(raw$results, c("Feature", "PC", "pvalue", "Significant"),
               ignore.order = TRUE)
  expect_named(adj$results, c("Feature", "PC", "pvalue", "pvalue_adj", "Significant"),
               ignore.order = TRUE)

  # pvalue always stays raw, regardless of whether p_adj was requested
  ok <- !is.na(raw$results$pvalue) & !is.na(adj$results$pvalue)
  expect_true(any(ok))
  expect_equal(adj$results$pvalue[ok], raw$results$pvalue[ok])

  # pvalue_adj is >= the raw p value
  ok_adj <- ok & !is.na(adj$results$pvalue_adj)
  expect_true(all(adj$results$pvalue_adj[ok_adj] >= adj$results$pvalue[ok_adj] - 1e-12))
})

test_that("ds.plotDrivers (split) includes pvalue_adj when p_adj is set", {

  setup <- setup_two_site_dslite()
  conns <- setup$conns
  withr::defer(DSI::datashield.logout(conns))

  out <- ds.plotDrivers(
    data = "expr", vars = "clinical", datasources = conns,
    type = "split", n_pc = 3L,
    p_adj = "BH", return_pvalues = TRUE
  )

  # Adjustment happens independently per site for type = "split"
  for (site_res in out$results) {
    expect_named(site_res,
                 c("Feature", "PC", "pvalue", "pvalue_adj", "Significant", "Site"),
                 ignore.order = TRUE)
    ok <- !is.na(site_res$pvalue) & !is.na(site_res$pvalue_adj)
    expect_true(all(site_res$pvalue_adj[ok] >= site_res$pvalue[ok] - 1e-12))
  }
})

test_that("ds.plotDrivers does not double-adjust p values with a single datasource", {

  setup <- setup_single_site_dslite()
  conns <- setup$conns
  withr::defer(DSI::datashield.logout(conns))

  # type = "combine" with a single site: combining is the identity operation,
  # and the adjustment happens once, server-side (see server_p_adj in
  # ds.plotDrivers). Applying p.adjust() again here on the already-adjusted
  # pvalue_adj must reproduce the same values (not a further-inflated result).
  out <- ds.plotDrivers(
    data = "expr", vars = "clinical", datasources = conns,
    type = "combine", n_pc = 3L,
    p_adj = "BH", return_pvalues = TRUE
  )

  res <- out$results
  expect_named(res, c("Feature", "PC", "pvalue", "pvalue_adj", "Significant"),
               ignore.order = TRUE)

  ok <- !is.na(res$pvalue)
  expected_adj <- stats::p.adjust(res$pvalue[ok], method = "BH")
  expect_equal(res$pvalue_adj[ok], expected_adj)
})

test_that("ds.plotDrivers (combine, n_pc = 5) reproduces known p values for site_data", {

  setup <- setup_two_site_dslite()
  conns <- setup$conns
  withr::defer(DSI::datashield.logout(conns))

  out <- ds.plotDrivers(
    data           = "expr",
    vars           = "clinical",
    datasources    = conns,
    type           = "combine",
    n_pc           = 5L,
    sig_cutoff     = 0.05,
    return_pvalues = TRUE
  )

  res <- out$results

  # 6 clinical variables x 5 PCs
  expect_equal(nrow(res), 30L)

  # No p_adj was requested, so there is no pvalue_adj column.
  expect_false("pvalue_adj" %in% colnames(res))

  get_row <- function(feature, pc) res[res$Feature == feature & res$PC == pc, ]

  # disease_status is an overwhelming driver of PC1 (p ~ 1e-29)
  row <- get_row("disease_status", "PC1")
  expect_lt(row$pvalue, 1e-20)
  expect_true(row$Significant)

  # age is not associated with PC1 (p ~ 0.26)
  row <- get_row("age", "PC1")
  expect_equal(row$pvalue, 0.26, tolerance = 0.02)
  expect_false(row$Significant)

  # batch is essentially unrelated to PC1 (p ~ 0.85)
  row <- get_row("batch", "PC1")
  expect_equal(row$pvalue, 0.85, tolerance = 0.05)
  expect_false(row$Significant)

  # age is borderline-significant on PC4 (p ~ 0.046, just under the 0.05 cutoff).
  # An absolute range check is used here rather than expect_equal(tolerance = ...):
  # testthat/waldo tolerance is relative, so on a value this small a tolerance
  # wide enough to be useful would barely constrain anything.
  row <- get_row("age", "PC4")
  expect_gt(row$pvalue, 0.03)
  expect_lt(row$pvalue, 0.05)
  expect_true(row$Significant)
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
