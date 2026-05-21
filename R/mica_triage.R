#' Classify MICA Cells by Reportability
#'
#' Applies simple reportability labels to MICA diagnostics so users can
#' distinguish tightly identified, prior-dominant, and point-estimate-only
#' cells.
#'
#' @param diagnostics A diagnostics data frame from `mica_diagnostics()`.
#' @param pd_pinned_cutoff Width cutoff for PD-pinned cells.
#' @param regression_pinned_cutoff Width cutoff for regression-pinned cells.
#' @param data_dominant_cutoff Minimum anchor `R^2` threshold for
#'   data-dominant reportability.
#' @param prior_dominant_cutoff Maximum anchor `R^2` threshold for
#'   prior-dominant reportability.
#' @param wide_pd_cutoff Minimum PD width used to label a cell as wide enough to
#'   be honestly prior-dominant.
#'
#' @return The input diagnostics data frame with `flag` and `reportability`
#'   columns added.
#' @export
mica_triage <- function(diagnostics,
                        pd_pinned_cutoff = 0.05,
                        regression_pinned_cutoff = 0.10,
                        data_dominant_cutoff = 0.7,
                        prior_dominant_cutoff = 0.3,
                        wide_pd_cutoff = 0.3) {
  if (!is.data.frame(diagnostics)) {
    stop("`diagnostics` must be a data frame returned by `mica_diagnostics()`.",
         call. = FALSE)
  }

  out <- diagnostics
  out$flag <- "observed"
  out$reportability <- NA_character_

  missing_rows <- which(!out$observed)
  for (row in missing_rows) {
    pd_pinned <- !is.na(out$width_pd[row]) && out$width_pd[row] < pd_pinned_cutoff
    reg_pinned <- !is.na(out$width_pred[row]) &&
      out$width_pred[row] < regression_pinned_cutoff
    out$flag[row] <- if (pd_pinned) {
      "PD-pinned"
    } else if (reg_pinned) {
      "regression-pinned"
    } else {
      "floating"
    }

    r2_min <- suppressWarnings(min(out$r2_i[row], out$r2_j[row], na.rm = TRUE))
    if (!is.finite(r2_min)) {
      out$reportability[row] <- "point_only"
    } else if (r2_min > data_dominant_cutoff) {
      out$reportability[row] <- "interval_data_dominant"
    } else if (r2_min < prior_dominant_cutoff &&
               (is.na(out$width_pd[row]) || out$width_pd[row] > wide_pd_cutoff)) {
      out$reportability[row] <- "interval_prior_dominant"
    } else {
      out$reportability[row] <- "point_only"
    }
  }

  out
}
