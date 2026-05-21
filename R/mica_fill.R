#' Complete a Partial Correlation Matrix with a PD-Projected Fill
#'
#' Applies a deterministic fill routine to a partial correlation matrix and
#' projects the result to a positive-definite correlation matrix.
#'
#' This is a practical deterministic MICA baseline. It is a
#' regression-based, PD-projected fill and should not be described as the
#' formal GJSW maximum-determinant completion.
#'
#' @param corr_matrix A symmetric correlation matrix with ones on the diagonal
#'   and `NA` entries for missing off-diagonal cells.
#' @param max_iter Maximum number of fill iterations.
#' @param tol Convergence tolerance on the maximum absolute cell update.
#' @param damping Damping factor applied to each update.
#' @param pd_tol Positive-definiteness tolerance used during projection.
#'
#' @return A list with the completed matrix and fill diagnostics.
#' @export
mica_fill_pd <- function(corr_matrix,
                         max_iter = 200,
                         tol = 1e-8,
                         damping = 0.5,
                         pd_tol = 1e-4) {
  validated <- .mica_validate_matrix(corr_matrix)
  x <- validated$corr_matrix
  missing_idx <- .mica_upper_missing(x)

  if (nrow(missing_idx) == 0) {
    completed <- x
    diag(completed) <- 1
    return(list(
      matrix = completed,
      iterations = 0L,
      converged = TRUE,
      max_delta = 0,
      n_missing = 0L,
      min_eigenvalue = min(eigen(completed, symmetric = TRUE,
                                 only.values = TRUE)$values)
    ))
  }

  work <- x
  work[is.na(work)] <- 0
  diag(work) <- 1
  work <- .mica_near_correlation(work, posd_tol = pd_tol)
  observed_mask <- !is.na(x)
  work[observed_mask] <- x[observed_mask]
  diag(work) <- 1

  iteration_done <- 0L
  delta_done <- Inf
  for (iteration in seq_len(max_iter)) {
    previous <- work
    for (row in seq_len(nrow(missing_idx))) {
      i <- missing_idx[row, "i"]
      j <- missing_idx[row, "j"]
      others <- setdiff(seq_len(nrow(work)), c(i, j))
      r_oo <- work[others, others, drop = FALSE]
      r_i <- work[i, others]
      r_j <- work[j, others]
      proposal <- as.numeric(t(r_i) %*% .mica_safe_solve(r_oo) %*% r_j)
      proposal <- .mica_clip_r(proposal)
      updated <- damping * proposal + (1 - damping) * work[i, j]
      work[i, j] <- updated
      work[j, i] <- updated
    }

    work <- .mica_near_correlation(work, posd_tol = pd_tol)
    work[observed_mask] <- x[observed_mask]
    diag(work) <- 1

    delta_done <- max(abs(work - previous), na.rm = TRUE)
    iteration_done <- iteration
    if (is.finite(delta_done) && delta_done < tol) {
      break
    }
  }

  eigenvalues <- eigen(work, symmetric = TRUE, only.values = TRUE)$values
  list(
    matrix = work,
    iterations = iteration_done,
    converged = is.finite(delta_done) && delta_done < tol,
    max_delta = delta_done,
    n_missing = nrow(missing_idx),
    min_eigenvalue = min(eigenvalues)
  )
}
