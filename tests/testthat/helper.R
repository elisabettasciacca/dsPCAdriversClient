# Test helpers for dsPCAdriversClient -----------------------------------------
#
# Shared utilities used across all test files. `testthat` sources this file
# automatically before running any test.
#
# The end-to-end tests use DSLite to spin up an in-memory two-site DataSHIELD
# environment, reusing the `site_data` example dataset shipped with the
# package. Tests degrade gracefully when DSLite (or any of its server-side
# companion packages) is not installed.

# Skip helpers ----------------------------------------------------------------

skip_if_no_dslite <- function() {
  testthat::skip_if_not_installed("DSLite")
  testthat::skip_if_not_installed("DSI")
  testthat::skip_if_not_installed("dsBase")
  testthat::skip_if_not_installed("dsBaseClient")
  testthat::skip_if_not_installed("dsSwissKnife")
  testthat::skip_if_not_installed("dsSwissKnifeClient")
  testthat::skip_if_not_installed("dsPCAdrivers")

  # dsSwissKnifeClient calls some DSI functions (e.g. datashield.aggregate)
  # without namespacing them, so DSI must be attached to the search path.
  suppressPackageStartupMessages(library(DSI))
}

# Mock DSConnection for unit tests --------------------------------------------

#' Build a fake DSConnection object that satisfies methods::is(., "DSConnection").
#' Used to exercise input-validation helpers without a real backend.
fake_ds_connection <- function() {
  testthat::skip_if_not_installed("DSI")
  # setClass is idempotent — calling it multiple times in a session is safe.
  if (!methods::isClass("FakeDSConn")) {
    methods::setClass("FakeDSConn", contains = "DSConnection")
  }
  methods::new("FakeDSConn")
}

# DSLite session --------------------------------------------------------------

#' Set permissive DSLite privacy options (NEVER use these with real data).
.set_dslite_test_options <- function() {
  options(
    datashield.privacyControlLevel = "permissive",
    nfilter.tab     = 3,
    nfilter.subset  = 3,
    nfilter.glm     = 1.0,
    nfilter.string  = 80,
    nfilter.kNN     = 3
  )
}

#' Configure a DSLite server with dsBase + dsSwissKnife methods and the
#' plotDriversDS aggregate method registered.
.configure_dslite_server <- function(dslite_server) {
  dslite_server$config(
    DSLite::defaultDSConfiguration(include = c("dsBase", "dsSwissKnife"))
  )
  dslite_server$aggregateMethod("plotDriversDS", dsPCAdrivers::plotDriversDS)
  dslite_server
}

#' Create a two-site DSLite session preloaded with site_data.
#'
#' Returns a list with the connection object and the loaded data. The caller
#' is responsible for `DSI::datashield.logout()` — use
#' `withr::defer(DSI::datashield.logout(conns))` inside each test.
#'
#' The DSLite server object is assigned to .GlobalEnv under the name
#' "dslite_server" because the DSLiteDriver looks the server up by name from
#' the calling frames. `withr::defer()` removes it after the test.
setup_two_site_dslite <- function() {
  skip_if_no_dslite()

  e <- new.env()
  utils::data("site_data", package = "dsPCAdriversClient", envir = e)
  site_data <- e$site_data

  dslite_server <- .configure_dslite_server(DSLite::newDSLiteServer(
    tables = list(
      site1_expr     = site_data$site1$expr,
      site2_expr     = site_data$site2$expr,
      site1_clinical = site_data$site1$clinical,
      site2_clinical = site_data$site2$clinical
    )
  ))

  assign("dslite_server", dslite_server, envir = globalenv())
  withr::defer(
    suppressWarnings(rm("dslite_server", envir = globalenv())),
    envir = parent.frame()
  )

  .set_dslite_test_options()

  builder <- DSI::newDSLoginBuilder()
  builder$append(server = "site1", url = "dslite_server",
                 table  = "site1_expr", driver = "DSLiteDriver")
  builder$append(server = "site2", url = "dslite_server",
                 table  = "site2_expr", driver = "DSLiteDriver")

  conns <- DSI::datashield.login(logins = builder$build(), assign = FALSE)

  DSI::datashield.assign.table(
    conns, "expr",
    c(site1 = "site1_expr", site2 = "site2_expr")
  )
  DSI::datashield.assign.table(
    conns, "clinical",
    c(site1 = "site1_clinical", site2 = "site2_clinical")
  )

  list(conns = conns, site_data = site_data)
}

#' Single-site variant of the DSLite session, useful for testing the
#' "single server forces type = split" branch of ds.plotDrivers().
setup_single_site_dslite <- function() {
  skip_if_no_dslite()

  e <- new.env()
  utils::data("site_data", package = "dsPCAdriversClient", envir = e)
  site_data <- e$site_data

  dslite_server <- .configure_dslite_server(DSLite::newDSLiteServer(
    tables = list(
      site1_expr     = site_data$site1$expr,
      site1_clinical = site_data$site1$clinical
    )
  ))

  assign("dslite_server", dslite_server, envir = globalenv())
  withr::defer(
    suppressWarnings(rm("dslite_server", envir = globalenv())),
    envir = parent.frame()
  )

  .set_dslite_test_options()

  builder <- DSI::newDSLoginBuilder()
  builder$append(server = "site1", url = "dslite_server",
                 table  = "site1_expr", driver = "DSLiteDriver")
  conns <- DSI::datashield.login(logins = builder$build(), assign = FALSE)

  DSI::datashield.assign.table(conns, "expr",     c(site1 = "site1_expr"))
  DSI::datashield.assign.table(conns, "clinical", c(site1 = "site1_clinical"))

  list(conns = conns, site_data = site_data)
}

# Synthetic site_results / site_metadata for combine_pvalues_fisher tests -----

#' Build minimal `site_results` + `site_metadata` structures matching what
#' plotDriversDS would return for two sites.
make_fake_site_outputs <- function(pvals_site1, pvals_site2,
                                   features = c("age", "sex"),
                                   pcs      = c("PC1", "PC2")) {

  build_one <- function(pvals) {
    combs <- expand.grid(Feature = features, PC = pcs,
                         stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    data.frame(
      Feature = combs$Feature,
      PC      = combs$PC,
      pvalue  = pvals,
      stringsAsFactors = FALSE
    )
  }

  site_results <- list(
    site1 = build_one(pvals_site1),
    site2 = build_one(pvals_site2)
  )
  site_metadata <- list(
    site1 = list(pc_names = pcs, var_names = features, n_observations = 60L),
    site2 = list(pc_names = pcs, var_names = features, n_observations = 50L)
  )

  list(site_results = site_results, site_metadata = site_metadata)
}
