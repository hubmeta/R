# Internal helpers for MICA input validation and matrix operations.

.mica_ridge <- 1e-8

.mica_clip_r <- function(r, eps = 1e-6) {
  pmax(pmin(r, 1 - eps), -1 + eps)
}

.mica_safe_solve <- function(x, ridge = .mica_ridge) {
  out <- try(solve(x), silent = TRUE)
  if (inherits(out, "try-error")) {
    out <- solve(x + diag(ridge, nrow(x)))
  }
  out
}

.mica_is_pd <- function(x, tol = 1e-8) {
  if (!isTRUE(all.equal(x, t(x), tolerance = 1e-6))) {
    return(FALSE)
  }
  values <- try(eigen(x, symmetric = TRUE, only.values = TRUE)$values,
                silent = TRUE)
  if (inherits(values, "try-error")) {
    return(FALSE)
  }
  all(values > tol)
}

.mica_near_correlation <- function(x, posd_tol = 1e-4, eps = 1e-6) {
  out <- as.matrix(Matrix::nearPD(
    x,
    corr = TRUE,
    keepDiag = TRUE,
    ensureSymmetry = TRUE,
    posd.tol = posd_tol
  )$mat)
  eigenvalues <- eigen(out, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalues) < posd_tol) {
    out <- out + diag(eps, nrow(out))
    out <- out / sqrt(outer(diag(out), diag(out)))
    diag(out) <- 1
  }
  out
}

.mica_upper_missing <- function(x) {
  idx <- which(is.na(x) & upper.tri(x), arr.ind = TRUE)
  colnames(idx) <- c("i", "j")
  idx
}

.mica_validate_matrix <- function(corr_matrix, k_matrix = NULL,
                                  tau_matrix = NULL) {
  x <- as.matrix(corr_matrix)
  storage.mode(x) <- "numeric"

  if (!is.numeric(x) || length(dim(x)) != 2 || nrow(x) != ncol(x)) {
    stop("`corr_matrix` must be a numeric square matrix.", call. = FALSE)
  }
  if (!isTRUE(all.equal(x, t(x), tolerance = 1e-8, check.attributes = FALSE))) {
    stop("`corr_matrix` must be symmetric.", call. = FALSE)
  }
  if (any(!is.na(diag(x)) & abs(diag(x) - 1) > 1e-8)) {
    stop("Diagonal entries of `corr_matrix` must equal 1.", call. = FALSE)
  }
  diag(x) <- 1
  if (any(!is.na(x) & abs(x) > 1)) {
    stop("Off-diagonal correlations must lie in [-1, 1].", call. = FALSE)
  }

  validate_side_matrix <- function(mat, name) {
    if (is.null(mat)) {
      return(NULL)
    }
    out <- as.matrix(mat)
    storage.mode(out) <- "numeric"
    if (!identical(dim(out), dim(x))) {
      stop(sprintf("`%s` must have the same dimensions as `corr_matrix`.", name),
           call. = FALSE)
    }
    if (!isTRUE(all.equal(out, t(out), tolerance = 1e-8,
                          check.attributes = FALSE))) {
      stop(sprintf("`%s` must be symmetric.", name), call. = FALSE)
    }
    out
  }

  list(
    corr_matrix = x,
    k_matrix = validate_side_matrix(k_matrix, "k_matrix"),
    tau_matrix = validate_side_matrix(tau_matrix, "tau_matrix")
  )
}

.mica_predict_cell <- function(corr_matrix, i, j, ridge = 0) {
  others <- setdiff(seq_len(nrow(corr_matrix)), c(i, j))
  obs_i <- others[!is.na(corr_matrix[i, others])]
  obs_j <- others[!is.na(corr_matrix[j, others])]
  anchors <- intersect(obs_i, obs_j)

  if (length(anchors) < 2) {
    return(list(rhat = NA_real_, r2_i = NA_real_, r2_j = NA_real_,
                anchors = anchors))
  }

  r_kk <- corr_matrix[anchors, anchors, drop = FALSE]
  if (anyNA(r_kk)) {
    keep <- apply(r_kk, 1, function(row) !anyNA(row))
    anchors <- anchors[keep]
    if (length(anchors) < 2) {
      return(list(rhat = NA_real_, r2_i = NA_real_, r2_j = NA_real_,
                  anchors = anchors))
    }
    r_kk <- corr_matrix[anchors, anchors, drop = FALSE]
  }

  r_kk_solve <- if (ridge > 0) r_kk + diag(ridge, nrow(r_kk)) else r_kk
  r_kk_inv <- .mica_safe_solve(r_kk_solve)
  r_i <- corr_matrix[i, anchors]
  r_j <- corr_matrix[j, anchors]
  beta_i <- as.vector(r_kk_inv %*% r_i)
  beta_j <- as.vector(r_kk_inv %*% r_j)
  r2_i <- as.numeric(t(beta_i) %*% r_kk %*% beta_i)
  r2_j <- as.numeric(t(beta_j) %*% r_kk %*% beta_j)
  rhat <- as.numeric(t(beta_i) %*% r_kk %*% beta_j)

  list(
    rhat = .mica_clip_r(rhat),
    r2_i = r2_i,
    r2_j = r2_j,
    anchors = anchors
  )
}

.mica_pd_interval <- function(filled_matrix, i, j, step = 0.01) {
  work <- filled_matrix
  lower <- NA_real_
  upper <- NA_real_

  for (value in seq(-0.99, 0.99, by = step)) {
    work[i, j] <- value
    work[j, i] <- value
    if (.mica_is_pd(work)) {
      lower <- value
      break
    }
  }
  for (value in seq(0.99, -0.99, by = -step)) {
    work[i, j] <- value
    work[j, i] <- value
    if (.mica_is_pd(work)) {
      upper <- value
      break
    }
  }

  c(lower = lower, upper = upper)
}
