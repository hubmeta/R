#' Matrix Imputation with Correlation-Aware Diagnostics
#'
#' Completes a partial correlation matrix with a deterministic
#' positive-definite fill and attaches per-cell diagnostics and triage labels.
#'
#' This is the package-facing MICA entry point for the current development
#' version of `hubmeta`. It returns a positive-definite completed matrix plus
#' reportability diagnostics. Bayesian uncertainty layers can be added on top
#' of this object in later versions.
#'
#' @param corr_matrix A symmetric correlation matrix with ones on the diagonal
#'   and `NA` entries for missing off-diagonal cells.
#' @param k_matrix Optional symmetric study-count matrix aligned to
#'   `corr_matrix`.
#' @param tau_matrix Optional symmetric heterogeneity matrix aligned to
#'   `corr_matrix`.
#' @param pd_step Step size for PD-feasibility interval sweeps.
#' @param ... Additional arguments passed to `mica_fill_pd()`.
#'
#' @return An object of class `"mica"`.
#' @export
mica <- function(corr_matrix,
                 k_matrix = NULL,
                 tau_matrix = NULL,
                 pd_step = 0.01,
                 ...) {
  validated <- .mica_validate_matrix(corr_matrix, k_matrix, tau_matrix)
  fill <- mica_fill_pd(validated$corr_matrix, ...)
  diagnostics <- mica_diagnostics(
    corr_matrix = validated$corr_matrix,
    k_matrix = validated$k_matrix,
    tau_matrix = validated$tau_matrix,
    fill = fill,
    pd_step = pd_step
  )
  triage <- mica_triage(diagnostics)

  structure(
    list(
      input_matrix = validated$corr_matrix,
      k_matrix = validated$k_matrix,
      tau_matrix = validated$tau_matrix,
      diagnostics = diagnostics,
      triage = triage,
      fill = fill,
      completed_matrix = fill$matrix,
      method = "deterministic_pd_fill",
      call = match.call()
    ),
    class = "mica"
  )
}
