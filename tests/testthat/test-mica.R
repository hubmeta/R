test_that("mica_fill_pd returns a positive-definite completion", {
  toy <- matrix(c(
    1.00, 0.32,   NA, 0.21,
    0.32, 1.00, 0.28,   NA,
      NA, 0.28, 1.00, 0.41,
    0.21,   NA, 0.41, 1.00
  ), nrow = 4, byrow = TRUE)

  fit <- mica_fill_pd(toy)

  expect_true(is.matrix(fit$matrix))
  expect_equal(diag(fit$matrix), rep(1, 4))
  expect_true(all(eigen(fit$matrix, symmetric = TRUE)$values > 0))
})

test_that("mica_diagnostics and triage produce expected columns", {
  toy <- matrix(c(
    1.00, 0.32,   NA, 0.21,
    0.32, 1.00, 0.28,   NA,
      NA, 0.28, 1.00, 0.41,
    0.21,   NA, 0.41, 1.00
  ), nrow = 4, byrow = TRUE)

  fit <- mica(toy)

  expect_s3_class(fit, "mica")
  expect_true(all(c("flag", "reportability") %in% names(fit$triage)))
  expect_true(any(!fit$triage$observed))
  expect_true(all(!is.na(fit$triage$flag[!fit$triage$observed])))
})

test_that("mica accepts the stan engine", {
  skip_on_cran()

  toy <- matrix(c(
    1.00, 0.32,   NA, 0.21,
    0.32, 1.00, 0.28,   NA,
      NA, 0.28, 1.00, 0.41,
    0.21,   NA, 0.41, 1.00
  ), nrow = 4, byrow = TRUE)

  fit <- mica(
    toy,
    engine = "stan",
    chains = 1,
    iter = 120,
    warmup = 60,
    refresh = 0,
    seed = 123
  )

  expect_s3_class(fit, "mica_stan")
  expect_true(is.matrix(as.matrix(fit)))
  expect_equal(dim(as.matrix(fit)), c(4, 4))
})
