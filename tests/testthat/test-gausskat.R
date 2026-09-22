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

  parallel_fit <- GausSKAT(
    Z,
    null_model,
    number_grid_points = 5L,
    warn_on_endpoint = FALSE,
    parallel = TRUE,
    n_cores = 2L
  )

  expect_true(parallel_fit$parallel)
  expect_equal(parallel_fit$n.cores, 2L)
  expect_equal(parallel_fit$ell.grid, fit$ell.grid, tolerance = 0)
  expect_equal(
    parallel_fit$component.p.values,
    fit$component.p.values,
    tolerance = 1e-12
  )
  expect_equal(parallel_fit$p.value, fit$p.value, tolerance = 1e-12)

  D <- weighted_squared_distances(Z, fit$weights)
  psock_fits <- GausSKAT:::.evaluate_component_fits(
    ell_grid = fit$ell.grid,
    Z = Z,
    null_model = null_model,
    weights = fit$weights,
    method = "davies",
    parallel_plan = list(workers = 2L, backend = "PSOCK")
  )
  psock_p_values <- vapply(
    psock_fits,
    function(component_fit) component_fit$p.value,
    numeric(1)
  )
  expect_equal(psock_p_values, fit$component.p.values, tolerance = 1e-12)

  reference_p_values <- vapply(fit$ell.grid, function(ell) {
    reference_kernel <- exp(-D / ell)
    reference_fit <- SKAT::SKAT(
      Z = Z,
      obj = null_model,
      kernel = reference_kernel,
      method = "davies",
      weights = rep(1, ncol(Z)),
      max_maf = 1
    )
    as.numeric(reference_fit$p.value)
  }, numeric(1))
  expect_equal(
    fit$component.p.values,
    reference_p_values,
    tolerance = 1e-10
  )
})

test_that("binary null models are rejected explicitly", {
  Z <- rbind(c(0, 1), c(1, 0), c(2, 1), c(0, 0))
  null_model <- structure(list(out_type = "D"), class = "SKAT_NULL_Model")
  expect_error(GausSKAT(Z, null_model), "continuous-trait")
})
