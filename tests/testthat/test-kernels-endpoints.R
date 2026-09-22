test_that("distance and kernel match the original weighted formula", {
  Z <- rbind(
    c(0, 0, 1),
    c(0, 1, 1),
    c(2, 0, 0),
    c(1, 2, 1)
  )
  weights <- c(25, 4, 1)
  weighted_Z <- sweep(Z, 2L, weights, `*`)
  reference_D <- as.matrix(stats::dist(weighted_Z))^2

  D <- weighted_squared_distances(Z, weights)
  expect_equal(D, reference_D, tolerance = 1e-12)
  expect_equal(
    gausskat_kernel(D, 100),
    exp(-reference_D / 100),
    tolerance = 1e-12
  )
})

test_that("lower endpoint applies the positive-distance median heuristic", {
  Z <- rbind(
    c(0, 0),
    c(0, 0),
    c(1, 0),
    c(2, 1)
  )
  D <- weighted_squared_distances(Z, c(2, 1))
  X <- cbind(1, c(-1, 0, 1, 2))
  grid <- make_adaptive_grid(D, X, number_grid_points = 5L)
  positive <- D[upper.tri(D) & D > 0]

  expect_equal(grid$ell_min, stats::median(positive))
  expect_gte(
    mean(exp(-positive / grid$ell_min) >= exp(-1) - 1e-14),
    0.5
  )
  expect_length(grid$ell_grid, 5L)
  expect_true(all(diff(grid$ell_grid) >= 0))
})

test_that("adaptive endpoints are equivariant to common effective-weight scale", {
  set.seed(2)
  Z <- matrix(sample(0:2, 60, replace = TRUE), nrow = 10)
  X <- cbind(1, stats::rnorm(10))
  weights <- seq(1, 2, length.out = ncol(Z))
  scale_factor <- 3

  first <- make_adaptive_grid(weighted_squared_distances(Z, weights), X)
  second <- make_adaptive_grid(
    weighted_squared_distances(Z, sqrt(scale_factor) * weights),
    X
  )

  expect_equal(second$ell_min, scale_factor * first$ell_min,
               tolerance = 1e-10)
  expect_equal(second$ell_max_approx,
               scale_factor * first$ell_max_approx,
               tolerance = 1e-9)
})

test_that("large hyperparameters approach the projected weighted-linear shape", {
  set.seed(3)
  Z <- matrix(sample(0:2, 80, replace = TRUE), nrow = 10)
  weights <- seq(1, 3, length.out = ncol(Z))
  D <- weighted_squared_distances(Z, weights)
  X <- cbind(1, stats::rnorm(nrow(Z)))
  projection <- GausSKAT:::.make_covariate_projection(X)
  gaussian_shape <- GausSKAT:::.scaled_projected_kernel(
    ell = max(D) * 1e7,
    D = D,
    projection = projection
  )
  weighted_Z <- sweep(Z, 2L, weights, `*`)
  linear_projected <- GausSKAT:::.project_symmetric_matrix(
    tcrossprod(weighted_Z),
    projection
  )
  linear_shape <- linear_projected / sqrt(sum(linear_projected^2))

  expect_lt(sqrt(sum((gaussian_shape - linear_shape)^2)), 1e-5)
})
