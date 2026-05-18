## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.width = 7,
  fig.height = 5
)

# DSLite demo chunks run only if all required packages are available
run_dslite <- requireNamespace("DSLite",            quietly = TRUE) &&
              requireNamespace("DSI",               quietly = TRUE) &&
              requireNamespace("dsBaseClient",       quietly = TRUE) &&
              requireNamespace("dsSwissKnifeClient", quietly = TRUE) &&
              requireNamespace("dsPCAdrivers",       quietly = TRUE)

## ----install, eval = FALSE----------------------------------------------------
# # Install the client package
# devtools::install_github("elisabettasciacca/dsPCAdriversClient")
# 
# # The server-side package must be installed on each DataSHIELD server
# # by the server administrator:
# # devtools::install_github("elisabettasciacca/dsPCAdrivers")

## ----dslite-setup, eval = run_dslite------------------------------------------
library(DSLite)
library(DSI)
library(dsPCAdrivers)
library(dsPCAdriversClient)
library(dsBaseClient)
library(dsSwissKnifeClient)
library(dsPCAdriversClient)

# Load simulated example data (2 sites: 60 and 50 samples, 30 genes each)
data(site_data)

# Create a DSLite server holding all four tables
dslite.server <- DSLite::newDSLiteServer(
  tables = list(
    site1_expr     = site_data$site1$expr,
    site2_expr     = site_data$site2$expr,
    site1_clinical = site_data$site1$clinical,
    site2_clinical = site_data$site2$clinical
  )
)

# Load dsBase and dsSwissKnife methods, then register the server-side function
dslite.server$config(DSLite::defaultDSConfiguration(include = c("dsBase", "dsSwissKnife")))
dslite.server$aggregateMethod("plotDriversDS", dsPCAdrivers::plotDriversDS)

# Relax privacy filters for local testing only.
# These settings must NEVER be used with real sensitive data.
# nfilter.glm = 1.0 is required here because the gene/sample ratio (30/60 = 0.5)
# exceeds the default threshold of 0.33 used by ds.cov internally.
options(
  datashield.privacyControlLevel = "permissive",
  nfilter.tab     = 3,
  nfilter.subset  = 3,
  nfilter.glm     = 1.0,
  nfilter.string  = 80,
  nfilter.kNN     = 3
)

## ----dslite-login, eval = run_dslite------------------------------------------
builder <- DSI::newDSLoginBuilder()
builder$append(server = "site1", url = "dslite.server",
               table = "site1_expr", driver = "DSLiteDriver")
builder$append(server = "site2", url = "dslite.server",
               table = "site2_expr", driver = "DSLiteDriver")

conns <- DSI::datashield.login(logins = builder$build(), assign = FALSE)

# Assign expression matrices (samples x genes — already in the correct orientation)
DSI::datashield.assign.table(
  conns  = conns,
  symbol = "expression_data",
  table  = c(site1 = "site1_expr", site2 = "site2_expr")
)

# Assign clinical metadata
DSI::datashield.assign.table(
  conns  = conns,
  symbol = "clinical_vars",
  table  = c(site1 = "site1_clinical", site2 = "site2_clinical")
)

