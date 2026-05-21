#' @export
print.mica <- function(x, ...) {
  missing_rows <- sum(!x$triage$observed)
  reportable <- sum(x$triage$reportability %in%
                      c("interval_data_dominant", "interval_prior_dominant"),
                    na.rm = TRUE)
  cat("<mica>\n")
  cat("  method: ", x$method, "\n", sep = "")
  cat("  variables: ", nrow(x$completed_matrix), "\n", sep = "")
  cat("  missing cells completed: ", missing_rows, "\n", sep = "")
  cat("  reportable interval cells: ", reportable, "\n", sep = "")
  invisible(x)
}

#' @export
print.mica_stan <- function(x, ...) {
  missing_rows <- nrow(x$posterior_cells)
  cat("<mica_stan>\n")
  cat("  variables: ", nrow(x$completed_matrix), "\n", sep = "")
  cat("  missing cells modeled: ", missing_rows, "\n", sep = "")
  cat("  chains: ", x$fit@sim$chains, "\n", sep = "")
  cat("  iter per chain: ", x$fit@sim$iter, "\n", sep = "")
  invisible(x)
}

#' @export
summary.mica <- function(object, ...) {
  triage <- object$triage
  structure(list(
    method = object$method,
    variables = nrow(object$completed_matrix),
    missing_cells = sum(!triage$observed),
    flags = table(triage$flag, useNA = "ifany"),
    reportability = table(triage$reportability, useNA = "ifany"),
    fill = object$fill
  ), class = "summary.mica")
}

#' @export
summary.mica_stan <- function(object, ...) {
  structure(list(
    variables = nrow(object$completed_matrix),
    missing_cells = nrow(object$posterior_cells),
    flags = table(object$triage$flag, useNA = "ifany"),
    reportability = table(object$triage$reportability, useNA = "ifany"),
    posterior_cells = object$posterior_cells
  ), class = "summary.mica_stan")
}

#' @export
print.summary.mica <- function(x, ...) {
  cat("MICA summary\n")
  cat("  method: ", x$method, "\n", sep = "")
  cat("  variables: ", x$variables, "\n", sep = "")
  cat("  missing cells: ", x$missing_cells, "\n", sep = "")
  cat("  fill iterations: ", x$fill$iterations, "\n", sep = "")
  cat("  min eigenvalue: ", signif(x$fill$min_eigenvalue, 4), "\n", sep = "")
  cat("\nFlags\n")
  print(x$flags)
  cat("\nReportability\n")
  print(x$reportability)
  invisible(x)
}

#' @export
print.summary.mica_stan <- function(x, ...) {
  cat("MICA Stan summary\n")
  cat("  variables: ", x$variables, "\n", sep = "")
  cat("  missing cells: ", x$missing_cells, "\n", sep = "")
  cat("\nFlags\n")
  print(x$flags)
  cat("\nReportability\n")
  print(x$reportability)
  invisible(x)
}

#' Extract the Completed Matrix from a MICA Object
#'
#' @param x A `"mica"` object.
#' @param ... Unused.
#'
#' @return A completed positive-definite correlation matrix.
#' @export
as.matrix.mica <- function(x, ...) {
  x$completed_matrix
}

#' @export
as.matrix.mica_stan <- function(x, ...) {
  x$completed_matrix
}
