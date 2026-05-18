# Generate example data for dsPCAdriversClient --------------------------------
#
# Produces a single list object `site_data` with two simulated sites,
# saved as data/site_data.rda. Kept intentionally small to comply with
# CRAN package size limits.
#
# Data structure:
#   site_data$site1$expr      data.frame: 60 samples x 30 genes (samples x genes)
#   site_data$site1$clinical  data.frame: 60 samples x 7 clinical variables
#   site_data$site2$expr      data.frame: 50 samples x 30 genes
#   site_data$site2$clinical  data.frame: 50 samples x 7 clinical variables
#
# Expression matrices are in samples x genes orientation, as required by
# dssPrincomp (which computes covariance across columns).
#
# Biological structure embedded in the data:
#   - Disease signal: cases upregulate the first 15 genes (+2 on log scale)
#   - Batch effect:   first 8 genes shifted by batch (batch 1/2/3 x 0.5)
#
# Example usage with DSLite:
#
#   data(site_data)
#
#   dslite.server <- DSLite::newDSLiteServer(
#     tables = list(
#       site1_expr     = site_data$site1$expr,
#       site2_expr     = site_data$site2$expr,
#       site1_clinical = site_data$site1$clinical,
#       site2_clinical = site_data$site2$clinical
#     )
#   )
#
# Run this script once (from the package root) to regenerate the .rda file.
# Do not source this file from package code.

set.seed(123)

n_features <- 30
site_sizes <- c(site1 = 60, site2 = 50)

# Helper: generate expression matrix and clinical data for one site ------------

create_site_data <- function(n, n_features) {

  expr_mat <- matrix(
    rnorm(n_features * n, mean = 10, sd = 2),
    nrow     = n_features,
    ncol     = n,
    dimnames = list(
      paste0("gene_",   seq_len(n_features)),
      paste0("sample_", seq_len(n))
    )
  )

  clinical_df <- data.frame(
    sample_id      = paste0("sample_", seq_len(n)),
    age            = rnorm(n, mean = 50, sd = 10),
    sex            = sample(c("M", "F"), n, replace = TRUE),
    batch          = factor(sample(1:3, n, replace = TRUE)),
    disease_status = sample(c("case", "control"), n, replace = TRUE),
    bmi            = rnorm(n, mean = 25, sd = 5),
    smoking        = sample(c("never", "former", "current"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )

  # Disease signal: cases upregulate the first 15 genes
  case_idx <- which(clinical_df$disease_status == "case")
  if (length(case_idx) > 0) {
    expr_mat[1:15, case_idx] <- expr_mat[1:15, case_idx] + 2
  }

  # Batch effect on the first 8 genes
  for (b in seq_len(3)) {
    batch_idx <- which(clinical_df$batch == b)
    if (length(batch_idx) > 0) {
      expr_mat[1:8, batch_idx] <- expr_mat[1:8, batch_idx] + b * 0.5
    }
  }

  list(
    expr     = as.data.frame(t(expr_mat)),  # samples x genes
    clinical = clinical_df
  )
}

# Generate and collect into a named list ---------------------------------------

site_data <- Map(create_site_data, n = site_sizes,
                 MoreArgs = list(n_features = n_features))

# Save to data/ ----------------------------------------------------------------

usethis::use_data(site_data, overwrite = TRUE)
