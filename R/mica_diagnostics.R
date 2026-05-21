#' Compute MICA Per-Cell Diagnostics
#'
#' Computes anchor-regression and PD-feasibility diagnostics for each
#' upper-triangle cell in a partial correlation matrix.
#'
#' @param corr_matrix A symmetric correlation matrix with ones on the diagonal
#'   and `NA` entries for missing off-diagonal cells.
#' @param k_matrix Optional matrix of study counts aligned to `corr_matrix`.
#' @param tau_matrix Optional matrix of heterogeneity values aligned to
#'   `corr_matrix`.
#' @param fill Optional output from `mica_fill_pd()`. If omitted, a
#'   deterministic fill is computed internally.
#' @param pd_step Step size used when sweeping the PD-feasibility interval.
#'
#' @return A data frame with one row per upper-triangle cell.
#' @export
mica_diagnostics <- function(corr_matrix,
                             k_matrix = NULL,
                             tau_matrix = NULL,
                             fill = NULL,
                             pd_step = 0.01) {
  validated <- .mica_validate_matrix(corr_matrix, k_matrix, tau_matrix)
  x <- validated$corr_matrix
  k_mat <- validated$k_matrix
  tau_mat <- validated$tau_matrix

  if (is.null(fill)) {
    fill <- mica_fill_pd(x)
  }

  filled_matrix <- fill$matrix
  idx <- which(upper.tri(x), arr.ind = TRUE)
  if (is.null(rownames(x))) {
    rownames(x) <- colnames(x) <- paste0("V", seq_len(nrow(x)))
  }

  diagnostics <- data.frame(
    i = idx[, 1],
    j = idx[, 2],
    var_i = rownames(x)[idx[, 1]],
    var_j = rownames(x)[idx[, 2]],
    observed = !is.na(x[idx]),
    r_obs = x[idx],
    k_obs = if (!is.null(k_mat)) k_mat[idx] else NA_real_,
    tau_obs = if (!is.null(tau_mat)) tau_mat[idx] else NA_real_,
    rhat = NA_real_,
    r2_i = NA_real_,
    r2_j = NA_real_,
    width_pred = NA_real_,
    width_pd_lower = NA_real_,
    width_pd_upper = NA_real_,
    width_pd = NA_real_,
    stringsAsFactors = FALSE
  )

  for (row in seq_len(nrow(diagnostics))) {
    i <- diagnostics$i[row]
    j <- diagnostics$j[row]
    pred <- .mica_predict_cell(x, i, j)

    diagnostics$rhat[row] <- pred$rhat
    diagnostics$r2_i[row] <- pred$r2_i
    diagnostics$r2_j[row] <- pred$r2_j
    if (!is.na(pred$r2_i) && !is.na(pred$r2_j)) {
      diagnostics$width_pred[row] <- sqrt(max(
        0,
        (1 - pred$r2_i) * (1 - pred$r2_j)
      ))
    }

    if (!diagnostics$observed[row]) {
      interval <- .mica_pd_interval(filled_matrix, i, j, step = pd_step)
      diagnostics$width_pd_lower[row] <- interval["lower"]
      diagnostics$width_pd_upper[row] <- interval["upper"]
      if (!is.na(interval["lower"]) && !is.na(interval["upper"])) {
        diagnostics$width_pd[row] <- interval["upper"] - interval["lower"]
      }
    } else {
      diagnostics$width_pd[row] <- 0
    }
  }

  diagnostics
}
