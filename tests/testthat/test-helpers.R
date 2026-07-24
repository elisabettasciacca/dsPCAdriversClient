# Unit tests for internal helper functions of dsPCAdriversClient -------------
#
# These tests exercise pure-R helpers and do NOT require DSLite or any
# DataSHIELD backend. Non-exported functions are accessed via `:::`.

# .check_args -----------------------------------------------------------------

test_that(".check_args rejects NULL, non-character and length != 1", {

  expect_error(
    dsPCAdriversClient:::.check_args(NULL, "data"),
    "Argument 'data' must not be NULL"
  )
  expect_error(
    dsPCAdriversClient:::.check_args(42, "data"),
    "must be a character string"
  )
  expect_error(
    dsPCAdriversClient:::.check_args(c("a", "b"), "data"),
    "must have length 1"
  )
})

test_that(".check_args accepts a single character string silently", {
  expect_silent(dsPCAdriversClient:::.check_args("expression", "data"))
})


# .check_datasources ----------------------------------------------------------

test_that(".check_datasources rejects non-list inputs and lists with wrong types", {

  expect_error(
    dsPCAdriversClient:::.check_datasources("not_a_list"),
    "must be a list of DSConnection objects"
  )
  expect_error(
    dsPCAdriversClient:::.check_datasources(list("string", 1)),
    "must be a list of DSConnection objects"
  )
})

test_that(".check_datasources accepts a list of DSConnection-class objects", {
  testthat::skip_if_not_installed("DSI")

  conn <- fake_ds_connection()
  expect_silent(dsPCAdriversClient:::.check_datasources(list(conn)))
  expect_silent(dsPCAdriversClient:::.check_datasources(list(site1 = conn, site2 = conn)))
})


# validate_federated_pcs ------------------------------------------------------

test_that("validate_federated_pcs passes when all sites share the same PC names", {

  meta <- list(
    site1 = list(pc_names = paste0("PC", 1:5)),
    site2 = list(pc_names = paste0("PC", 1:5)),
    site3 = list(pc_names = paste0("PC", 1:5))
  )
  expect_true(dsPCAdriversClient:::validate_federated_pcs(meta))
})

test_that("validate_federated_pcs errors when PC names differ across sites", {

  meta <- list(
    site1 = list(pc_names = paste0("PC", 1:5)),
    site2 = list(pc_names = paste0("PC", 1:4))
  )
  expect_error(
    dsPCAdriversClient:::validate_federated_pcs(meta),
    "PC names differ across sites"
  )
})


# combine_pvalues_fisher ------------------------------------------------------

test_that("combine_pvalues_fisher matches the closed-form Fisher formula", {

  # Standard Fisher: chi^2 = -2 * sum(log(p)), df = 2k
  p1 <- 0.01
  p2 <- 0.04
  expected_chi  <- -2 * (log(p1) + log(p2))
  expected_pval <- stats::pchisq(expected_chi, df = 4, lower.tail = FALSE)

  fakes <- make_fake_site_outputs(
    pvals_site1 = c(p1, 0.5, 0.5, 0.5),
    pvals_site2 = c(p2, 0.5, 0.5, 0.5)
  )

  combined <- dsPCAdriversClient:::combine_pvalues_fisher(
    fakes$site_results, fakes$site_metadata
  )

  age_pc1 <- combined$pvalue[combined$Feature == "age" & combined$PC == "PC1"]
  expect_equal(age_pc1, expected_pval, tolerance = 1e-10)
})

test_that("combine_pvalues_fisher returns the single available p when others are NA", {

  fakes <- make_fake_site_outputs(
    pvals_site1 = c(0.02, 0.5, 0.5, 0.5),
    pvals_site2 = c(NA,   0.5, 0.5, 0.5)
  )
  combined <- dsPCAdriversClient:::combine_pvalues_fisher(
    fakes$site_results, fakes$site_metadata
  )
  age_pc1 <- combined$pvalue[combined$Feature == "age" & combined$PC == "PC1"]
  expect_equal(age_pc1, 0.02)
})

test_that("combine_pvalues_fisher returns NA when all sites are NA", {

  fakes <- make_fake_site_outputs(
    pvals_site1 = c(NA, 0.5, 0.5, 0.5),
    pvals_site2 = c(NA, 0.5, 0.5, 0.5)
  )
  combined <- dsPCAdriversClient:::combine_pvalues_fisher(
    fakes$site_results, fakes$site_metadata
  )
  age_pc1 <- combined$pvalue[combined$Feature == "age" & combined$PC == "PC1"]
  expect_true(is.na(age_pc1))
})

test_that("combine_pvalues_fisher produces one row per Feature x PC combination", {

  fakes <- make_fake_site_outputs(
    pvals_site1 = c(0.1, 0.2, 0.3, 0.4),
    pvals_site2 = c(0.1, 0.2, 0.3, 0.4)
  )
  combined <- dsPCAdriversClient:::combine_pvalues_fisher(
    fakes$site_results, fakes$site_metadata
  )
  expect_equal(nrow(combined), 4)
  expect_named(combined, c("Feature", "PC", "pvalue"), ignore.order = TRUE)
})


# format_pval_sci --------------------------------------------------------------

test_that("format_pval_sci formats p values as scientific notation with no zero-padded exponent", {
  expect_equal(dsPCAdriversClient:::format_pval_sci(3.123e-4), "3.12e-4")
  expect_equal(dsPCAdriversClient:::format_pval_sci(9.892336e-29), "9.89e-29")
  expect_equal(dsPCAdriversClient:::format_pval_sci(0.05), "5.00e-2")
  expect_true(is.na(dsPCAdriversClient:::format_pval_sci(NA_real_)))
})

# build_drivers_plot ----------------------------------------------------------

#' Minimal results frame matching what ds.plotDrivers produces internally.
make_results_df <- function(with_site = FALSE) {
  res <- expand.grid(
    Feature = c("age", "sex"),
    PC      = paste0("PC", 1:3),
    stringsAsFactors = FALSE
  )
  set.seed(1)
  res$pvalue      <- stats::runif(nrow(res))
  res$pvalue_adj  <- res$pvalue / 2  # arbitrary but distinguishable from pvalue
  res$Significant <- res$pvalue <= 0.5
  if (with_site) res$Site <- rep(c("site1", "site2"),
                                 length.out = nrow(res))
  res
}

test_that("build_drivers_plot returns a ggplot object", {
  testthat::skip_if_not_installed("ggplot2")

  res <- make_results_df()
  p <- dsPCAdriversClient:::build_drivers_plot(
    results        = res,
    sig_cutoff     = 0.05,
    p_adj          = NULL,
    max_col        = NULL,
    title          = "Test plot",
    legend         = "right",
    transpose_plot = FALSE,
    label          = FALSE,
    faceted        = FALSE
  )

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Test plot")
})

test_that("build_drivers_plot honours transpose_plot", {
  testthat::skip_if_not_installed("ggplot2")

  res <- make_results_df()

  p_normal <- dsPCAdriversClient:::build_drivers_plot(
    res, sig_cutoff = 0.05, p_adj = NULL, max_col = NULL,
    title = "", legend = "right",
    transpose_plot = FALSE, label = FALSE
  )
  p_transposed <- dsPCAdriversClient:::build_drivers_plot(
    res, sig_cutoff = 0.05, p_adj = NULL, max_col = NULL,
    title = "", legend = "right",
    transpose_plot = TRUE, label = FALSE
  )

  # In the non-transposed plot, x = PC and y = Feature.
  # In the transposed plot, the assignment is swapped.
  expect_equal(as.character(p_normal$data$x),     as.character(res$PC))
  expect_equal(as.character(p_transposed$data$x), as.character(res$Feature))
})

test_that("build_drivers_plot adds facets when faceted = TRUE", {
  testthat::skip_if_not_installed("ggplot2")

  res <- make_results_df(with_site = TRUE)
  p <- dsPCAdriversClient:::build_drivers_plot(
    results        = res,
    sig_cutoff     = 0.05,
    p_adj          = NULL,
    max_col        = NULL,
    title          = "Faceted plot",
    legend         = "right",
    transpose_plot = FALSE,
    label          = FALSE,
    faceted        = TRUE
  )
  expect_s3_class(p, "ggplot")
  # The facet attribute should be a FacetWrap (or a subclass of it)
  expect_true(inherits(p$facet, "FacetWrap"))
})

test_that("build_drivers_plot adds a geom_text layer when label = TRUE", {
  testthat::skip_if_not_installed("ggplot2")

  res <- make_results_df()

  p_no_label <- dsPCAdriversClient:::build_drivers_plot(
    res, sig_cutoff = 0.05, p_adj = NULL, max_col = NULL,
    title = "", legend = "right",
    transpose_plot = FALSE, label = FALSE
  )
  p_label <- dsPCAdriversClient:::build_drivers_plot(
    res, sig_cutoff = 0.05, p_adj = NULL, max_col = NULL,
    title = "", legend = "right",
    transpose_plot = FALSE, label = TRUE
  )

  expect_equal(length(p_label$layers), length(p_no_label$layers) + 1L)
})

test_that("build_drivers_plot labels tiles with the formatted p value", {
  testthat::skip_if_not_installed("ggplot2")

  res <- make_results_df()
  p <- dsPCAdriversClient:::build_drivers_plot(
    res, sig_cutoff = 0.05, p_adj = NULL, max_col = NULL,
    title = "", legend = "right",
    transpose_plot = FALSE, label = TRUE
  )

  label_layer_data <- ggplot2::layer_data(p, length(p$layers))
  expect_equal(label_layer_data$label,
               dsPCAdriversClient:::format_pval_sci(res$pvalue))
})

test_that("build_drivers_plot uses p adj legend label when p_adj is set", {
  testthat::skip_if_not_installed("ggplot2")

  res <- make_results_df()

  p_raw <- dsPCAdriversClient:::build_drivers_plot(
    res, sig_cutoff = 0.05, p_adj = NULL, max_col = NULL,
    title = "", legend = "right",
    transpose_plot = FALSE, label = FALSE
  )
  p_adj <- dsPCAdriversClient:::build_drivers_plot(
    res, sig_cutoff = 0.05, p_adj = "BH", max_col = NULL,
    title = "", legend = "right",
    transpose_plot = FALSE, label = FALSE
  )

  # Fill scale name differs between adjusted and raw plots
  fill_raw <- p_raw$scales$get_scales("fill")$name
  fill_adj <- p_adj$scales$get_scales("fill")$name
  expect_false(identical(fill_raw, fill_adj))
})
