.skat_matrix_kernel <- function(Z, null_model, K, method) {
  fit <- SKAT::SKAT(
    Z = Z,
    obj = null_model,
    kernel = K,
    method = method,
    weights = rep(1, ncol(Z)),
    max_maf = 1
  )

  if (length(fit$p.value) != 1L || !is.finite(fit$p.value)) {
    stop("SKAT did not return a finite component p-value.", call. = FALSE)
  }
  fit
}

#' GausSKAT test for a continuous trait
#'
#' Constructs a phenotype-independent, data-adaptive grid of weighted Gaussian
#' kernels, evaluates each kernel using the official `SKAT` package, and
#' aggregates the dependent component p-values by ACAT.
#'
#' The default per-variant weights are beta-density weights with parameters
#' `(1, 25)`, as in weighted SKAT. Supplied weights follow the `SKAT` package
#' convention: the effective coefficient in both the squared weighted distance
#' and the weighted-linear kernel is the square of the supplied weight.
#'
#' @param Z Complete genotype dosage matrix, with individuals in rows and
#'   variants in columns.
#' @param null_model A continuous-trait null model returned by
#'   [SKAT::SKAT_Null_Model()].
#' @param X Optional null-model design matrix. By default, `null_model$X1` is
#'   used. If supplied, it must describe the same retained observations and
#'   contain an intercept in its column space.
#' @param weights Optional non-negative per-variant weights.
#' @param weights_beta Beta-weight parameters used when `weights` is `NULL`.
#' @param epsilon Stabilization tolerance for the practical upper-endpoint
#'   approximation.
#' @param number_grid_points Number of logarithmically spaced kernels.
#' @param method SKAT p-value method passed to [SKAT::SKAT()].
#' @param acat_weights Optional ACAT component weights. Equal weights are used
#'   by default.
#' @param keep_component_fits Whether to retain the complete SKAT fit for every
#'   grid point.
#' @param warn_on_endpoint Whether to warn when the observed doubling change at
#'   the approximate upper endpoint exceeds `epsilon`.
#' @return An object of class `GausSKAT` containing the aggregated p-value,
#'   component p-values, adaptive grid, and endpoint diagnostics.
#' @export
GausSKAT <- function(Z, null_model, X = NULL, weights = NULL,
                     weights_beta = c(1, 25), epsilon = 0.05,
                     number_grid_points = 5L, method = "davies",
                     acat_weights = NULL, keep_component_fits = FALSE,
                     warn_on_endpoint = TRUE) {
  call <- match.call()
  Z <- .validate_genotypes(Z)

  if (!inherits(null_model, "SKAT_NULL_Model")) {
    stop("null_model must be returned by SKAT::SKAT_Null_Model().",
         call. = FALSE)
  }
  if (is.null(null_model$out_type) || null_model$out_type != "C") {
    stop("GausSKAT currently supports continuous-trait null models only.",
         call. = FALSE)
  }
  if (!method %in% c("davies", "liu", "liu.mod")) {
    stop("method must be 'davies', 'liu', or 'liu.mod'.", call. = FALSE)
  }

  included_rows <- .included_rows(Z, null_model)
  Z_analysis <- Z[included_rows, , drop = FALSE]
  X_analysis <- .resolve_design_matrix(
    X = X,
    null_model = null_model,
    included_rows = included_rows,
    n_total = nrow(Z)
  )

  maf <- .minor_allele_frequencies(Z_analysis)
  if (is.null(weights)) {
    weights <- .beta_weights(maf, weights_beta)
  } else {
    weights <- as.numeric(weights)
    if (length(weights) != ncol(Z)) {
      stop("weights must have one value for each column of Z.",
           call. = FALSE)
    }
    if (any(!is.finite(weights)) || any(weights < 0)) {
      stop("weights must be finite and non-negative.", call. = FALSE)
    }
  }

  D <- weighted_squared_distances(Z_analysis, weights)
  grid <- make_adaptive_grid(
    D = D,
    X = X_analysis,
    epsilon = epsilon,
    number_grid_points = number_grid_points
  )

  if (isTRUE(warn_on_endpoint) &&
      grid$delta_at_ell_max > epsilon * (1 + 1e-8)) {
    warning(
      paste0(
        "The first-order upper-endpoint approximation has an observed ",
        "doubling change of ", signif(grid$delta_at_ell_max, 4),
        ", which exceeds epsilon = ", signif(epsilon, 4), "."
      ),
      call. = FALSE
    )
  }

  component_fits <- lapply(grid$ell_grid, function(ell) {
    .skat_matrix_kernel(
      Z = Z,
      null_model = null_model,
      K = gausskat_kernel(D, ell),
      method = method
    )
  })
  component_p_values <- vapply(
    component_fits,
    function(fit) as.numeric(fit$p.value),
    numeric(1)
  )

  result <- list(
    p.value = acat_combine(component_p_values, acat_weights),
    component.p.values = component_p_values,
    ell.grid = grid$ell_grid,
    ell.min = grid$ell_min,
    ell.max.approx = grid$ell_max_approx,
    ell.max.first.order = grid$ell_max_first_order,
    delta.at.ell.max = grid$delta_at_ell_max,
    alignment.at.ell.max = grid$alignment_at_ell_max,
    epsilon = epsilon,
    grid.degenerate = grid$grid_degenerate,
    weights = weights,
    maf = maf,
    method = method,
    number.grid.points = number_grid_points,
    rank.X = grid$rank_X,
    call = call
  )

  if (isTRUE(keep_component_fits)) {
    result$component.fits <- component_fits
  }

  class(result) <- "GausSKAT"
  result
}

#' @export
print.GausSKAT <- function(x, ...) {
  cat("GausSKAT\n")
  cat("  ACAT p-value:       ", format.pval(x$p.value, digits = 6), "\n",
      sep = "")
  cat("  Grid points:        ", length(x$ell.grid), "\n", sep = "")
  cat("  ell_min:            ", format(x$ell.min, digits = 6), "\n",
      sep = "")
  cat("  ell_max,approx:     ", format(x$ell.max.approx, digits = 6),
      "\n", sep = "")
  cat("  Change at ell_max:  ", format(x$delta.at.ell.max, digits = 5),
      "\n", sep = "")
  invisible(x)
}
