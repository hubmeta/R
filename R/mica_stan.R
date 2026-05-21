 .mica_fit_stan <- function(corr_matrix,
                           k_matrix = NULL,
                           tau_matrix = NULL,
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
  validated <- .mica_validate_matrix(corr_matrix, k_matrix, tau_matrix)
  x <- validated$corr_matrix
  k_mat <- validated$k_matrix

  fill <- mica_fill_pd(x)
  diagnostics <- mica_diagnostics(
    corr_matrix = x,
    k_matrix = k_mat,
    tau_matrix = validated$tau_matrix,
    fill = fill
  )
  triage <- mica_triage(diagnostics)

  miss_idx <- .mica_upper_missing(x)
  obs_idx <- which(!is.na(x) & upper.tri(x), arr.ind = TRUE)
  colnames(obs_idx) <- c("i", "j")

  n_miss <- nrow(miss_idx)
  n_obs <- nrow(obs_idx)

  prior_mean_z <- numeric(n_miss)
  prior_prec_z <- numeric(n_miss)
  for (m in seq_len(n_miss)) {
    i <- miss_idx[m, "i"]
    j <- miss_idx[m, "j"]
    pred <- .mica_predict_cell(x, i, j)
    rhat <- if (is.na(pred$rhat)) fill$matrix[i, j] else pred$rhat
    rhat_clipped <- .mica_clip_r(rhat)
    if (!is.na(pred$r2_i) && !is.na(pred$r2_j)) {
      width_pred <- sqrt(max(0, (1 - pred$r2_i) * (1 - pred$r2_j)))
      sd_z <- width_pred / (1 - rhat_clipped^2)
      prior_prec_z[m] <- 1 / max(sd_z^2, 1e-3)
    } else {
      prior_prec_z[m] <- 1
    }
    prior_prec_z[m] <- max(prior_prec_floor, prior_prec_z[m] * prior_prec_mult)
    prior_mean_z[m] <- atanh(rhat_clipped)
  }

  z_obs <- numeric(n_obs)
  within_prec_obs <- numeric(n_obs)
  for (o in seq_len(n_obs)) {
    i <- obs_idx[o, "i"]
    j <- obs_idx[o, "j"]
    z_obs[o] <- atanh(.mica_clip_r(x[i, j]))
    k_val <- if (!is.null(k_mat)) k_mat[i, j] else NA_real_
    if (is.na(k_val) || k_val < 1) {
      k_val <- 1
    }
    within_prec_obs[o] <- max(1, k_val * (n_avg_per_study - 3))
  }

  standata <- list(
    d = nrow(x),
    n_obs = n_obs,
    obs_idx = unname(obs_idx),
    z_obs = z_obs,
    within_prec_obs = within_prec_obs,
    n_miss = n_miss,
    miss_idx = unname(miss_idx),
    prior_mean_z = prior_mean_z,
    prior_prec_z = prior_prec_z,
    tau2_scale = tau2_scale,
    lkj_eta = lkj_eta
  )

  fit <- rstan::sampling(
    stanmodels$mica_corr,
    data = standata,
    chains = chains,
    iter = iter,
    warmup = warmup,
    seed = seed,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    refresh = refresh,
    ...
  )

  summary_fit <- rstan::summary(fit, pars = "R_out")$summary
  completed_mean <- x
  completed_q025 <- x
  completed_q975 <- x

  cell_summary <- data.frame(
    i = integer(0),
    j = integer(0),
    mean = numeric(0),
    sd = numeric(0),
    q025 = numeric(0),
    q975 = numeric(0)
  )

  if (n_miss > 0) {
    for (m in seq_len(n_miss)) {
      i <- miss_idx[m, "i"]
      j <- miss_idx[m, "j"]
      row_name <- sprintf("R_out[%d,%d]", i, j)
      vals <- summary_fit[row_name, c("mean", "sd", "2.5%", "97.5%")]
      completed_mean[i, j] <- completed_mean[j, i] <- vals["mean"]
      completed_q025[i, j] <- completed_q025[j, i] <- vals["2.5%"]
      completed_q975[i, j] <- completed_q975[j, i] <- vals["97.5%"]
      cell_summary <- rbind(
        cell_summary,
        data.frame(
          i = i,
          j = j,
          mean = vals["mean"],
          sd = vals["sd"],
          q025 = vals["2.5%"],
          q975 = vals["97.5%"]
        )
      )
    }
  }

  structure(
    list(
      input_matrix = x,
      k_matrix = k_mat,
      tau_matrix = validated$tau_matrix,
      deterministic_fill = fill,
      diagnostics = diagnostics,
      triage = triage,
      standata = standata,
      fit = fit,
      posterior_cells = cell_summary,
      completed_matrix = completed_mean,
      completed_q025 = completed_q025,
      completed_q975 = completed_q975,
      call = match.call()
    ),
    class = "mica_stan"
  )
}

#' Bayesian MICA Completion with Stan
#'
#' Compatibility wrapper for users who want to explicitly request the Stan
#' engine. Internally this forwards to `mica(..., engine = "stan")`.
#'
#' @inheritParams mica
#' @param chains Number of Stan chains.
#' @param iter Number of post-warmup iterations per chain.
#' @param warmup Number of warmup iterations per chain.
#' @param seed Random seed passed to Stan.
#' @param adapt_delta Stan `adapt_delta`.
#' @param max_treedepth Stan `max_treedepth`.
#' @param prior_prec_floor Minimum prior precision on the Fisher-z scale.
#' @param prior_prec_mult Multiplier applied to the regression-bound prior
#'   precision.
#' @param tau2_scale Scale parameter for the half-Cauchy prior on `tau`.
#' @param lkj_eta LKJ concentration parameter.
#' @param n_avg_per_study Average per-study sample size used when only `k`
#'   counts are available.
#' @param refresh Stan refresh interval.
#' @param ... Additional arguments passed to [rstan::sampling()].
#'
#' @return An object of class `"mica_stan"`.
#' @export
#'
#' @examples
#' \donttest{
#' toy <- matrix(c(
#'   1.00, 0.32,   NA, 0.21,
#'   0.32, 1.00, 0.28,   NA,
#'     NA, 0.28, 1.00, 0.41,
#'   0.21,   NA, 0.41, 1.00
#' ), nrow = 4, byrow = TRUE)
#'
#' fit <- mica_stan(toy, chains = 1, iter = 200, warmup = 100, refresh = 0)
#' fit
#' }
mica_stan <- function(corr_matrix,
                      k_matrix = NULL,
                      tau_matrix = NULL,
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
  mica(
    corr_matrix = corr_matrix,
    k_matrix = k_matrix,
    tau_matrix = tau_matrix,
    engine = "stan",
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
  )
}
