#' Weighted squared genetic distances
#'
#' Computes pairwise squared Euclidean distances after multiplying genotype
#' column `k` by `weights[k]`. Accordingly, the effective coefficient of that
#' variant in the squared distance is `weights[k]^2`, matching the weight
#' convention used by weighted-linear SKAT.
#'
#' @param Z Complete genotype dosage matrix with individuals in rows.
#' @param weights Non-negative per-variant weights supplied to SKAT.
#' @return A symmetric matrix of weighted squared genetic distances.
#' @export
weighted_squared_distances <- function(Z, weights) {
  Z <- .validate_genotypes(Z)
  weights <- as.numeric(weights)

  if (length(weights) != ncol(Z)) {
    stop("weights must have one value for each column of Z.", call. = FALSE)
  }
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("weights must be finite and non-negative.", call. = FALSE)
  }
  if (!any(weights > 0)) {
    stop("At least one variant weight must be positive.", call. = FALSE)
  }

  weighted_Z <- sweep(Z, 2L, weights, `*`)
  squared_norm <- rowSums(weighted_Z * weighted_Z)
  D <- outer(squared_norm, squared_norm, `+`) - 2 * tcrossprod(weighted_Z)
  .validate_distance_matrix(D)
}

#' Construct a weighted Gaussian kernel
#'
#' @param D Matrix of weighted squared genetic distances.
#' @param ell Positive Gaussian-kernel hyperparameter.
#' @return The matrix with entries `exp(-D[ii']/ell)`.
#' @export
gausskat_kernel <- function(D, ell) {
  D <- .validate_distance_matrix(D)
  if (length(ell) != 1L || !is.finite(ell) || ell <= 0) {
    stop("ell must be a finite positive number.", call. = FALSE)
  }
  exp(-D / ell)
}

.ell_min_positive_median <- function(D) {
  positive_distances <- D[upper.tri(D) & D > 0]
  if (length(positive_distances) == 0L) {
    stop("No strictly positive off-diagonal genetic distances were found.",
         call. = FALSE)
  }
  stats::median(positive_distances)
}

.ell_max_first_order <- function(D, projection, epsilon) {
  A <- -.project_symmetric_matrix(D, projection)
  B <- 0.5 * .project_symmetric_matrix(D * D, projection)
  norm_A_squared <- sum(A * A)
  norm_B_squared <- sum(B * B)

  if (!is.finite(norm_A_squared) ||
      norm_A_squared <= 100 * .Machine$double.eps * max(1, norm_B_squared)) {
    stop("The first-order projected kernel has numerically zero norm.",
         call. = FALSE)
  }

  inner_AB <- sum(A * B)
  norm_B_perp_squared_raw <-
    norm_B_squared - inner_AB^2 / norm_A_squared
  roundoff_tolerance <-
    1000 * .Machine$double.eps * max(1, norm_B_squared)

  if (!is.finite(norm_B_perp_squared_raw) ||
      norm_B_perp_squared_raw < -roundoff_tolerance) {
    stop("The orthogonal second-order component is numerically invalid.",
         call. = FALSE)
  }

  norm_B_perp_squared <- max(norm_B_perp_squared_raw, 0)
  list(
    ell_max_first_order = sqrt(norm_B_perp_squared) /
      (2 * epsilon * sqrt(norm_A_squared)),
    norm_A = sqrt(norm_A_squared),
    norm_B = sqrt(norm_B_squared),
    norm_B_perp = sqrt(norm_B_perp_squared)
  )
}

.scaled_projected_kernel <- function(ell, D, projection) {
  if (length(ell) != 1L || !is.finite(ell) || ell <= 0) {
    stop("ell must be a finite positive number.", call. = FALSE)
  }

  # Since S 11' S = 0, expm1(-D/ell) has the same projection as
  # exp(-D/ell) and is more stable for very large ell.
  projected <- .project_symmetric_matrix(expm1(-D / ell), projection)
  projected_norm <- sqrt(sum(projected * projected))

  if (!is.finite(projected_norm) ||
      projected_norm <= .Machine$double.eps) {
    stop("The projected kernel has numerically zero Frobenius norm.",
         call. = FALSE)
  }
  projected / projected_norm
}

.kernel_change_at_ell <- function(ell, D, projection) {
  scaled_ell <- .scaled_projected_kernel(ell, D, projection)
  scaled_double <- .scaled_projected_kernel(2 * ell, D, projection)
  c(
    delta = sqrt(sum((scaled_double - scaled_ell)^2)),
    alignment = sum(scaled_ell * scaled_double)
  )
}

#' Construct the data-adaptive GausSKAT hyperparameter grid
#'
#' The lower endpoint is the median strictly positive weighted genetic
#' distance. The practical upper endpoint is the first-order approximation to
#' the point at which doubling the hyperparameter changes the Frobenius-scaled,
#' covariate-adjusted kernel by no more than `epsilon`.
#'
#' @param D Matrix of weighted squared genetic distances.
#' @param X Null-model design matrix. Its column space must contain an
#'   intercept.
#' @param epsilon Positive stabilization tolerance smaller than `sqrt(2)`.
#' @param number_grid_points Number of logarithmically spaced grid points.
#' @return A list containing the endpoints, grid, and diagnostics.
#' @export
make_adaptive_grid <- function(D, X, epsilon = 0.05,
                               number_grid_points = 5L) {
  D <- .validate_distance_matrix(D)
  X <- .validate_numeric_matrix(X, "X", minimum_rows = 2L)

  if (nrow(X) != nrow(D)) {
    stop("D and X must describe the same individuals.", call. = FALSE)
  }
  if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0 ||
      epsilon >= sqrt(2)) {
    stop("epsilon must be positive and smaller than sqrt(2).",
         call. = FALSE)
  }

  number_grid_points <- as.integer(number_grid_points)
  if (length(number_grid_points) != 1L || is.na(number_grid_points) ||
      number_grid_points < 2L) {
    stop("number_grid_points must be an integer of at least two.",
         call. = FALSE)
  }

  projection <- .make_covariate_projection(X)
  ell_min <- .ell_min_positive_median(D)
  upper <- .ell_max_first_order(D, projection, epsilon)
  grid_degenerate <- upper$ell_max_first_order <= ell_min
  ell_max_approx <- max(ell_min, upper$ell_max_first_order)
  ell_grid <- exp(seq(log(ell_min), log(ell_max_approx),
                      length.out = number_grid_points))
  diagnostic <- .kernel_change_at_ell(ell_max_approx, D, projection)

  list(
    ell_min = ell_min,
    ell_max_approx = ell_max_approx,
    ell_max_first_order = upper$ell_max_first_order,
    ell_grid = ell_grid,
    delta_at_ell_max = unname(diagnostic[["delta"]]),
    alignment_at_ell_max = unname(diagnostic[["alignment"]]),
    epsilon = epsilon,
    grid_degenerate = grid_degenerate,
    rank_X = projection$rank,
    number_columns_X = projection$number_columns,
    rank_deficient_X = projection$rank_deficient,
    intercept_error_X = projection$intercept_error
  )
}
