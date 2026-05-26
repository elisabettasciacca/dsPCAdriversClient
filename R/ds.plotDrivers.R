#' @title Plot drivers of variation in omic data
#'
#' @description Client-side function to visualise associations between principal
#' components and variables in a federated setting. PCA is computed internally
#' via \code{dsSwissKnifeClient::dssPrincomp} — no pre-existing PCA object is
#' required. Only aggregate statistics (p values) leave the servers.
#'
#' @param data Character string. Name of the data matrix or data frame on the
#'   server(s) on which PCA will be computed. The matrix must have
#'   \strong{samples as rows and features as columns}.
#'   Note that omic matrices (e.g. RNAseq)
#'   are typically stored in the opposite orientation (features x samples) and
#'   must be transposed before use (e.g. via \code{datashield.assign(conns,
#'   "expr_t", quote(t(expr)))}). Only numeric columns are used by
#'   \code{dssPrincomp}.
#' @param vars Character string. Name of the variables (predictors) data frame on the
#'   server(s) containing the features to associate with the PCs (e.g. clinical
#'   variables).
#' @param datasources A list of \code{DSConnection}-class objects. Default NULL
#'   uses all available connections.
#' @param type Character. Either \code{"combine"} or \code{"split"}.
#'   \itemize{
#'     \item \code{"combine"}: A federated (pooled) PCA is computed across all
#'       sites via \code{dssPrincomp(..., type = "combine")}. All sites share the
#'       same PC loadings. P values are then combined across sites using standard
#'       (unweighted) Fisher's method, yielding a single plot.
#'     \item \code{"split"}: A separate, site-specific PCA is computed on each
#'       server via \code{dssPrincomp(..., type = "split")}. Associations are
#'       computed independently per site and a faceted plot (one facet per site)
#'       is returned. P values are \emph{not} combined across sites.
#'   }
#'   Default \code{"combine"}.
#' @param center Logical. Should columns be centred before PCA? Default TRUE.
#' @param scale Logical. Should columns be scaled to unit variance before PCA?
#'   Default FALSE.
#' @param parametric Logical. Use parametric tests (Pearson correlation, ANOVA)?
#'   If FALSE, uses rank-based alternatives (Spearman, Kruskal-Wallis). Default
#'   TRUE.
#' @param n_pc Integer. Number of principal components to include. Default 5.
#' @param label Logical. Print association values on tiles? Default FALSE.
#' @param sig_cutoff Numeric. Significance threshold for outlining tiles. Default
#'   0.05.
#' @param p_adj Optional character. P value adjustment method: one of
#'   \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"},
#'   \code{"BH"}, \code{"BY"}, \code{"fdr"}. Default NULL (no adjustment). When
#'   \code{type = "combine"}, adjustment is applied after Fisher combination.
#' @param max_col Numeric. Maximum value for the colour scale. If NULL, uses the
#'   maximum observed \eqn{-\log_{10}(p)}. Default NULL.
#' @param title Character. Plot title. Default \code{"Variation By Feature"}.
#' @param legend Character. Legend position: one of \code{"bottom"},
#'   \code{"left"}, \code{"top"}, \code{"right"}, \code{"bottomright"},
#'   \code{"bottomleft"}, \code{"topleft"}, \code{"topright"}. Default
#'   \code{"right"}.
#' @param na_drop_threshold Integer. Minimum number of non-NA values required
#'   for a variable to be included. Default 4.
#' @param transpose_plot Logical. Transpose axes (PCs on y, features on x)?
#'   Default FALSE.
#' @param drop_insignificant_x Logical. Remove PCs with no significant
#'   associations? Default FALSE.
#' @param drop_insignificant_y Logical. Remove variables with no significant
#'   associations? Default FALSE.
#' @param return_pvalues Logical. Return results data frame(s) instead of a plot?
#'   Default FALSE.
#' @param verbose Logical. Print diagnostic messages from the server? Default
#'   FALSE.
#'
#' @return
#' When \code{return_pvalues = FALSE} (default):
#' \itemize{
#'   \item \code{type = "combine"}: a single \code{ggplot2} heatmap.
#'   \item \code{type = "split"}: a single \code{ggplot2} heatmap faceted by
#'     site.
#' }
#' When \code{return_pvalues = TRUE}, a list containing:
#' \itemize{
#'   \item \code{results}: data frame (or named list of data frames for
#'     \code{type = "split"} with multiple sites) with columns
#'     \code{Feature}, \code{PC}, \code{pvalue}, \code{Association},
#'     \code{Significant}.
#'   \item \code{metadata}: list with summary information (n per site, method,
#'     etc.).
#' }
#'
#' @details
#' \strong{PCA scores on the server}: After this function runs, each server will
#' retain an object named \code{paste0(data, "_scores")} (e.g. if
#' \code{data = "expr"}, the scores will be at \code{"expr_scores"}). This
#' object can be used for further downstream analyses (e.g. biplots via
#' \code{biplot.dssPrincomp}).
#'
#' \strong{Single datasource}: PCA is computed on that server and a single plot
#' is returned, regardless of \code{type}.
#'
#' \strong{type = "combine"}: \code{dssPrincomp} derives loadings from the
#' pooled covariance matrix (n-weighted across sites) and sends the same
#' loadings to every server. Each site projects its samples onto the shared PC
#' axes, making cross-site p value combination statistically meaningful. P values
#' are combined using standard Fisher's method (equal weight per site). Weighting
#' by sample size is deliberately avoided: because the pooled covariance matrix
#' is already n-weighted, each site's p values implicitly reflect its sample
#' size (a larger site produces smaller p values for the same effect size).
#' Applying an additional n-weight in the combination step would therefore
#' double-count the contribution of larger sites.
#'
#' \strong{type = "split"}: Each site computes its own PCA independently. The
#' resulting PCs may differ in direction and magnitude across sites, so p values
#' are \emph{not} combined. Results are shown as a faceted plot.
#'
#' Association strength is visualised as \eqn{-\log_{10}(p)}. Tiles are outlined
#' in black when p (or adjusted p) \eqn{\leq} \code{sig_cutoff}.
#'
#' @examples
#' \dontrun{
#' library(DSI)
#' library(DSOpal)
#' library(dsSwissKnifeClient)
#'
#' builder <- DSI::newDSLoginBuilder()
#' builder$append(server = "site1", url = "https://site1.org", ...)
#' builder$append(server = "site2", url = "https://site2.org", ...)
#' conns <- DSI::datashield.login(logins = builder$build(), assign = TRUE)
#'
#' # Example 1: Federated PCA (pooled) + combined p values ------------------
#' p <- ds.plotDrivers(
#'   data       = "expression",
#'   vars       = "clinical",
#'   datasources = conns,
#'   type       = "combine",
#'   center     = TRUE,
#'   scale      = TRUE,
#'   n_pc       = 10,
#'   sig_cutoff = 0.05
#' )
#' print(p)
#' # PC scores are now available on each server as "expression_scores"
#'
#' # Example 2: Site-specific PCA + separate plots per site -----------------
#' p <- ds.plotDrivers(
#'   data       = "expression",
#'   vars       = "clinical",
#'   datasources = conns,
#'   type       = "split",
#'   n_pc       = 5,
#'   p_adj      = "BH"
#' )
#' print(p)
#'
#' # Example 3: Return raw data instead of plot ------------------------------
#' res <- ds.plotDrivers(
#'   data        = "expression",
#'   vars        = "clinical",
#'   datasources = conns,
#'   type        = "combine",
#'   return_pvalues = TRUE
#' )
#' head(res$results)
#' }
#'
#' @export
ds.plotDrivers <- function(data = NULL,
                           vars = NULL,
                           datasources = NULL,
                           type = "combine",
                           center = TRUE,
                           scale = FALSE,
                           parametric = TRUE,
                           n_pc = 5L,
                           label = FALSE,
                           sig_cutoff = 0.05,
                           p_adj = NULL,
                           max_col = NULL,
                           title = "Variation By Feature",
                           legend = "right",
                           na_drop_threshold = 4,
                           transpose_plot = FALSE,
                           drop_insignificant_x = FALSE,
                           drop_insignificant_y = FALSE,
                           return_pvalues = FALSE,
                           verbose = FALSE) {

  # Validate input data
  if (is.null(datasources)) datasources <- DSI::datashield.connections_find()
  .check_datasources(datasources)

  .check_args(data, "data")
  .check_args(vars, "vars")

  if (!type %in% c("combine", "split")) {
    stop("'type' must be either 'combine' or 'split'.", call. = FALSE)
  }

  # Check that data and vars exist on all servers
  invisible(lapply(names(datasources), function(ds_name) {
    if (!check_object_exists(data, datasources[ds_name])) {
      stop("Object '", data, "' not found on server '", ds_name, "'.", call. = FALSE)
    }
    if (!check_object_exists(vars, datasources[ds_name])) {
      stop("Object '", vars, "' not found on server '", ds_name, "'.", call. = FALSE)
    }
  }))

  # Compute PCA via dssPrincomp
  # The scores object assigned on each server will be named paste0(data, "_scores").
  # This object persists after the function returns and can be used for further
  # downstream analyses
  scores_name <- paste0(data, "_scores")

  # With a single server, "type=combine" makes no practical difference,
  # use "split" to avoid unnecessary cross-server communication.
  pca_type <- if (length(datasources) == 1) "split" else type

  if (verbose) message("Running dssPrincomp (type = '", pca_type, "') ...")

  dsSwissKnifeClient::dssPrincomp(
    df          = data,
    type        = pca_type,
    center      = center,
    scale       = scale,
    datasources = datasources
  )

  if (verbose) message("PCA scores assigned to '", scores_name, "' on each server.")

  # Call server-side association function
  # When type = "combine", p value combination will be done here after aggregation,
  # (client-side) so we do not apply p_adj server-side.
  server_p_adj <- if (length(datasources) > 1 && type == "combine") NULL else p_adj

  cally <- call("plotDriversDS",
                pcs.name          = scores_name,
                vars.name         = vars,
                parametric        = parametric,
                n_pc              = as.integer(n_pc),
                na_drop_threshold = as.integer(na_drop_threshold),
                p_adj             = server_p_adj,
                verbose           = verbose)

  server_outputs <- DSI::datashield.aggregate(datasources, cally)

  # Parse server outputs
  site_results <- lapply(server_outputs, function(x) x$results)
  site_metadata <- lapply(server_outputs, function(x) {
    list(
      pc_names      = x$pc_names,
      var_names     = x$var_names,
      n_observations = x$n_observations
    )
  })


  if (type == "split") {

    if (length(datasources) > 1) {
      message("type = 'split': returning separate plots per site.")
    }

    site_results_labelled <- lapply(names(site_results), function(site_name) {
      res  <- site_results[[site_name]]
      meta <- site_metadata[[site_name]]

      res$Significant <- !is.na(res$pvalue) & res$pvalue <= sig_cutoff
      res$Feature     <- factor(res$Feature, levels = meta$var_names)
      res$PC          <- factor(res$PC,      levels = meta$pc_names)
      if (length(datasources) > 1) res$Site <- site_name
      return(res)
    })

    if (return_pvalues) {
      if (length(datasources) == 1) {
        return(list(
          results  = site_results_labelled[[1]],
          metadata = site_metadata[[1]]
        ))
      }
      return(list(
        results  = setNames(site_results_labelled, names(site_results)),
        metadata = site_metadata
      ))
    }

    if (length(datasources) == 1) {
      results <- site_results_labelled[[1]]
      if (drop_insignificant_x) {
        results <- results[results$PC %in% unique(results$PC[results$Significant]), ]
      }
      if (drop_insignificant_y) {
        results <- results[results$Feature %in% unique(results$Feature[results$Significant]), ]
      }
      return(build_drivers_plot(results, sig_cutoff, p_adj, max_col, title,
                                legend, transpose_plot, label))
    }

    combined_df <- do.call(rbind, site_results_labelled)

    if (drop_insignificant_x) {
      combined_df <- combined_df[
        combined_df$PC %in% unique(combined_df$PC[combined_df$Significant]), ]
    }
    if (drop_insignificant_y) {
      combined_df <- combined_df[
        combined_df$Feature %in% unique(combined_df$Feature[combined_df$Significant]), ]
    }

    return(build_drivers_plot(combined_df, sig_cutoff, p_adj, max_col,
                              title, legend, transpose_plot, label, faceted = TRUE))
  }

  # Multiple servers, type = "combine"
  if (length(datasources) > 1) {
    message("Combining associations across ", length(datasources), " sites ",
            "(Fisher's method).")
  }

  # Sanity check: with federated PCA, PC names must be identical across sites
  validate_federated_pcs(site_metadata)

  combined_results <- combine_pvalues_fisher(site_results, site_metadata)

  if (!is.null(p_adj)) {
    combined_results$pvalue      <- p.adjust(combined_results$pvalue, method = p_adj)
    combined_results$Association <- -log10(combined_results$pvalue)
  }

  combined_results$Significant <- !is.na(combined_results$pvalue) & combined_results$pvalue <= sig_cutoff
  combined_results$Feature     <- factor(combined_results$Feature,
                                         levels = site_metadata[[1]]$var_names)
  combined_results$PC          <- factor(combined_results$PC,
                                         levels = site_metadata[[1]]$pc_names)

  if (return_pvalues) {
    total_n <- sum(vapply(site_metadata, function(x) x$n_observations, integer(1)))
    return(list(
      results  = combined_results,
      metadata = list(
        n_sites        = length(datasources),
        total_n        = total_n,
        meta_method    = "fisher"
      )
    ))
  }

  if (drop_insignificant_x) {
    combined_results <- combined_results[
      combined_results$PC %in% unique(combined_results$PC[combined_results$Significant]), ]
  }
  if (drop_insignificant_y) {
    combined_results <- combined_results[
      combined_results$Feature %in% unique(combined_results$Feature[combined_results$Significant]), ]
  }

  if (title == "Variation By Feature") {
    total_n <- sum(vapply(site_metadata, function(x) x$n_observations, integer(1)))
    title   <- paste0("Variation By Feature (",
                      length(datasources), " sites, N=", total_n, ")")
  }

  return(build_drivers_plot(combined_results, sig_cutoff, p_adj, max_col, title,
                            legend, transpose_plot, label))
}

# Helper functions ---------------------------------------------------------------

.check_datasources <- function(datasources) {
  if (!(is.list(datasources) && all(unlist(lapply(datasources, function(d) {
    methods::is(d, "DSConnection")
  }))))) {
    stop("'datasources' must be a list of DSConnection objects.", call. = FALSE)
  }
}

.check_args <- function(x, arg_name) {
  if (is.null(x)) {
    stop("Argument '", arg_name, "' must not be NULL.", call. = FALSE)
  }
  if (!is.character(x)) {
    stop("Argument '", arg_name, "' must be a character string, not ",
         class(x), ".", call. = FALSE)
  }
  if (length(x) != 1) {
    stop("Argument '", arg_name, "' must have length 1.", call. = FALSE)
  }
}

check_object_exists <- function(obj_name, datasource) {
  symbols <- try(DSI::datashield.symbols(datasource), silent = TRUE)
  if (inherits(symbols, "try-error")) return(FALSE)
  # datashield.symbols() returns a named list (one element per server);
  # each element is a character vector of symbol names
  return(any(vapply(symbols, function(x) obj_name %in% x, logical(1))))
}

validate_federated_pcs <- function(site_metadata) {
  pc_names_list <- lapply(site_metadata, function(x) x$pc_names)
  reference_pcs <- pc_names_list[[1]]

  invisible(lapply(seq_along(pc_names_list)[-1], function(i) {
    if (!identical(pc_names_list[[i]], reference_pcs)) {
      stop("PC names differ across sites despite type = 'combine'. ",
           "This should not happen — please report this as a bug.", call. = FALSE)
    }
  }))
  return(invisible(TRUE))
}

combine_pvalues_fisher <- function(site_results, site_metadata) {

  all_pcs  <- site_metadata[[1]]$pc_names
  all_vars <- site_metadata[[1]]$var_names

  combinations <- expand.grid(
    Feature = all_vars,
    PC      = all_pcs,
    stringsAsFactors = FALSE
  )

  combinations$pvalue <- mapply(function(feat, pc) {

    pvals <- vapply(site_results, function(site_data) {
      idx <- which(site_data$Feature == feat & site_data$PC == pc)
      if (length(idx) == 0) return(NA_real_)
      site_data$pvalue[idx]
    }, numeric(1))

    valid_idx <- !is.na(pvals)
    pvals     <- pvals[valid_idx]

    if (length(pvals) == 0) return(NA_real_)
    if (length(pvals) == 1) return(pvals)

    # Standard (unweighted) Fisher's method. Additional n-weighting is deliberately
    # omitted: with type = "combine", the pooled covariance matrix used by
    # dssPrincomp is already n-weighted, so each site's p values implicitly
    # reflect its sample size. Weighting again here would double-count the
    # contribution of larger sites.
    fisher_stat <- -2 * sum(log(pvals))
    combined_p  <- pchisq(fisher_stat, df = 2 * length(pvals), lower.tail = FALSE)

    return(combined_p)

  }, feat = combinations$Feature, pc = combinations$PC)

  combinations$Association <- -log10(combinations$pvalue)
  return(combinations)
}


#' @title Build the PCA drivers heatmap
#'
#' @description Internal helper that constructs the ggplot heatmap of
#'   variable-PC associations used by \code{ds.plotDrivers}. Handles
#'   significance highlighting, p value adjustment annotations, axis
#'   transposition and optional faceting by study.
#'
#' @param results Data frame of associations returned by the aggregation step.
#' @param sig_cutoff Numeric. Significance threshold on -log10(p).
#' @param p_adj Character. Method used for p value adjustment (for labelling).
#' @param max_col Numeric or \code{NULL}. Upper limit for the colour scale.
#' @param title Character. Plot title.
#' @param legend Logical. Whether to display the legend.
#' @param transpose_plot Logical. If \code{TRUE}, swap features and PCs on the axes.
#' @param label Character. Label shown next to significant cells.
#' @param faceted Logical. If \code{TRUE}, facet the plot by study.
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_text coord_equal facet_wrap
#'  scale_fill_gradientn scale_colour_manual guides guide_legend labs theme_bw
#'  theme element_text element_rect
#' @keywords internal
build_drivers_plot <- function(results, sig_cutoff, p_adj, max_col,
                               title, legend, transpose_plot, label,
                               faceted = FALSE) {
  if (is.null(max_col)) {
    max_col <- max(ceiling(results$Association), na.rm = TRUE)
  }

  results$x <- if (transpose_plot) results$Feature else results$PC
  results$y <- if (transpose_plot) results$PC      else results$Feature

  leg_lab <- if (!is.null(p_adj)) {
    expression(-log[10](p[adj]))
  } else {
    expression(-log[10](p))
  }

  p <- ggplot(
    results,
    aes(x = x, y = y, fill = Association, colour = Significant)
  ) +
    geom_tile(linewidth = if (faceted) 0.5 else 1, width = 0.9, height = 0.9) +
    coord_equal() +
    scale_fill_gradientn(
      colours = c("white", "dodgerblue1", "dodgerblue3", "dodgerblue4"),
      name    = leg_lab,
      limits  = c(0, max_col)
    ) +
    scale_colour_manual(
      values = c("grey90", "black"),
      labels = c(
        paste(ifelse(is.null(p_adj), "p", "p adj"), ">",  sig_cutoff),
        paste(ifelse(is.null(p_adj), "p", "p adj"), "<=", sig_cutoff)
      ),
      name = ""
    ) +
    guides(
      colour = guide_legend(override.aes = list(fill = "white"))
    ) +
    labs(title = title, x = "", y = "") +
    theme_bw() +
    theme(
      plot.title      = element_text(hjust = 0.5),
      axis.text.x     = element_text(
        angle = ifelse(transpose_plot, 315, 0),
        hjust = ifelse(transpose_plot, 0, 0.5)),
      axis.text.y     = element_text(size = 11),
      legend.position = legend
    )

  if (faceted) {
    p <- p +
      facet_wrap(~ Site) +
      theme(
        strip.background = element_rect(fill = "#dce8f0", colour = "grey70"),
        strip.text       = element_text(face = "bold", size = 10),
        panel.spacing    = unit(1, "lines")
      )
  }

  if (label) {
    p <- p + geom_text(
      aes(label = round(Association, 2)),
      colour = "black",
      size   = if (faceted) 2.5 else 3
    )
  }

  return(p)
}
