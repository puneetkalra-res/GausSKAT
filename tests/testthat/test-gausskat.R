test_that("GausSKAT returns a complete five-component result", {
  skip_if_not_installed("SKAT")
  set.seed(11)
  n <- 60
  m <- 12
  Z <- matrix(stats::rbinom(n * m, size = 2, prob = 0.05), nrow = n)
  X1 <- stats::rnorm(n)
  X2 <- stats::rbinom(n, 1, 0.5)
  Y <- 0.5 * X1 + 0.5 * X2 + stats::rnorm(n)
  dat <- data.frame(Y = Y, X1 = X1, X2 = X2)
  null_model <- SKAT::SKAT_Null_Model(Y ~ X1 + X2, data = dat,
                                      out_type = "C")

  fit <- GausSKAT(
    Z,
    null_model,
    number_grid_points = 5L,
    warn_on_endpoint = FALSE
  )

  expect_s3_class(fit, "GausSKAT")
  expect_length(fit$ell.grid, 5L)
  expect_length(fit$component.p.values, 5L)
  expect_true(all(is.finite(fit$component.p.values)))
  expect_gte(fit$p.value, 0)
  expect_lte(fit$p.value, 1)
})

test_that("binary null models are rejected explicitly", {
  Z <- rbind(c(0, 1), c(1, 0), c(2, 1), c(0, 0))
  null_model <- structure(list(out_type = "D"), class = "SKAT_NULL_Model")
  expect_error(GausSKAT(Z, null_model), "continuous-trait")
})
