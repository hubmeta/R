#' Matrix Imputation with Correlation-Aware Diagnostics
#'
#' Completes a partial correlation matrix with either a deterministic
#' positive-definite fill or the Bayesian Stan engine, then attaches per-cell
#' diagnostics and triage labels.
#'
#' This is the package-facing MICA front door in `hubmeta`.
#'
#' @param corr_matrix A symmetric correlation matrix with ones on the diagonal
#'   and `NA` entries for missing off-diagonal cells.
#' @param k_matrix Optional symmetric study-count matrix aligned to
#'   `corr_matrix`.
#' @param tau_matrix Optional symmetric heterogeneity matrix aligned to
#'   `corr_matrix`.
#' @param engine Completion engine. Use `"deterministic"` for the
#'   PD-projected fill or `"stan"` for the Bayesian MICA model.
#' @param pd_step Step size for PD-feasibility interval sweeps.
#' @param chains Number of Stan chains when `engine = "stan"`.
#' @param iter Number of post-warmup iterations per chain when
#'   `engine = "stan"`.
#' @param warmup Number of warmup iterations per chain when `engine = "stan"`.
#' @param seed Random seed passed to Stan when `engine = "stan"`.
#' @param adapt_delta Stan `adapt_delta` when `engine = "stan"`.
#' @param max_treedepth Stan `max_treedepth` when `engine = "stan"`.
#' @param prior_prec_floor Minimum prior precision on the Fisher-z scale when
#'   `engine = "stan"`.
#' @param prior_prec_mult Multiplier applied to the regression-bound prior
#'   precision when `engine = "stan"`.
#' @param tau2_scale Scale parameter for the half-Cauchy prior on `tau` when
#'   `engine = "stan"`.
#' @param lkj_eta LKJ concentration parameter when `engine = "stan"`.
#' @param n_avg_per_study Average per-study sample size used when only `k`
#'   counts are available and `engine = "stan"`.
#' @param refresh Stan refresh interval when `engine = "stan"`.
#' @param ... Additional arguments passed to `mica_fill_pd()` or
#'   [rstan::sampling()] depending on `engine`.
#'
#' @return An object of class `"mica"` for deterministic fits or
#'   `"mica_stan"` for Bayesian fits.
#' @export
#'
#' @examples
#' toy <- matrix(c(
#'   1.00, 0.32,   NA, 0.21,
#'   0.32, 1.00, 0.28,   NA,
#'     NA, 0.28, 1.00, 0.41,
#'   0.21,   NA, 0.41, 1.00
#' ), nrow = 4, byrow = TRUE)
#'
#' det_fit <- mica(toy)
#' det_fit
#'
#' \donttest{
#' bayes_fit <- mica(toy, engine = "stan", chains = 1, iter = 200,
#'                   warmup = 100, refresh = 0)
#' bayes_fit
#' }
mica <- function(corr_matrix,
                 k_matrix = NULL,
                 tau_matrix = NULL,
                 engine = c("deterministic", "stan"),
                 pd_step = 0.01,
                 chains = 4,
                 iter = 1000,
                 warmup = floor(iter / 2),
                 seed = 20260521,
                 adapt_delta = 0.95,
                 max_treedepth = 12,
                 prior_prec_floor = 0.1,
                 prior_prec_mult = 1,
                 tau2_scale = 0.05,
                 lkj_eta = 1,
                 n_avg_per_study = 200,
                 refresh = 0,
                 ...) {
  engine <- match.arg(engine)
  validated <- .mica_validate_matrix(corr_matrix, k_matrix, tau_matrix)

  if (identical(engine, "stan")) {
    return(.mica_fit_stan(
      corr_matrix = validated$corr_matrix,
      k_matrix = validated$k_matrix,
      tau_matrix = validated$tau_matrix,
      chains = chains,
      iter = iter,
      warmup = warmup,
      seed = seed,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      prior_prec_floor = prior_prec_floor,
      prior_prec_mult = prior_prec_mult,
      tau2_scale = tau2_scale,
      lkj_eta = lkj_eta,
      n_avg_per_study = n_avg_per_study,
      refresh = refresh,
      ...
    ))
  }

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
      engine = engine,
      call = match.call()
    ),
    class = "mica"
  )
}
