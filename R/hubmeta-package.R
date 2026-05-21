#' @keywords internal
"_PACKAGE"

#' hubmeta: Meta-Analysis and Matrix Completion Tools
#'
#' Tools for psychometric meta-analysis and partial meta-analytic
#' correlation-matrix completion. The package currently includes
#' Hunter-Schmidt style pooling, Morris-weight meta-analysis, and MICA
#' utilities for diagnosing, deterministically completing, and Bayesianly
#' completing incomplete correlation matrices before downstream
#' multivariate modeling.
#'
#' @name hubmeta-package
#' @aliases hubmeta
#' @useDynLib hubmeta, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
NULL
