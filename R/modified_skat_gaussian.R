# This file isolates the continuous-trait weighted-Gaussian pathway used in
# the research implementation. The calculation is adapted from
# SKAT.linear.Other in SKAT and from the weighted-Gaussian modification in the
# author's research source. SKAT and this package are distributed under GPL
# (>= 2). Unrelated SKAT source is deliberately not copied into this package.

.modified_weighted_gaussian_kernel <- function(Z, ell, weights) {
  if (length(ell) != 1L || !is.finite(ell) || ell <= 0) {
    stop("ell must be a finite positive number.", call. = FALSE)
  }
  if (length(weights) != ncol(Z) || any(!is.finite(weights)) ||
      any(weights < 0)) {
    stop("weights must be finite, non-negative, and match ncol(Z).",
         call. = FALSE)
  }

  # This is the kernel.gaussian.weighted1 calculation used in the modified
  # SKAT source: genotype column k is multiplied by the supplied SKAT weight
  # before squared Euclidean distances are calculated.
  weighted_Z <- sweep(Z, 2L, weights, `*`)
  distance_squared <- as.matrix(stats::dist(
    weighted_Z,
    method = "euclidean"
  ))^2
  exp(-distance_squared / ell)
}

.skat_mixture_pvalue <- function(Q, W, Q_resampling, method) {
  function_name <- switch(
    method,
    davies = "Get_Davies_PVal",
    liu = "Get_Liu_PVal",
    liu.mod = "Get_Liu_PVal.MOD",
    stop("Unsupported SKAT p-value method.", call. = FALSE)
  )
  pvalue_function <- getFromNamespace(function_name, "SKAT")
  pvalue_function(Q, W, Q_resampling)
}

.modified_skat_gaussian <- function(Z, null_model, ell, weights,
                                    method = "davies") {
  if (!inherits(null_model, "SKAT_NULL_Model") ||
      is.null(null_model$out_type) || null_model$out_type != "C") {
    stop("The modified weighted-Gaussian pathway requires a continuous-trait SKAT null model.",
         call. = FALSE)
  }
  if (!method %in% c("davies", "liu", "liu.mod")) {
    stop("method must be 'davies', 'liu', or 'liu.mod'.", call. = FALSE)
  }

  Z <- .validate_genotypes(Z)
  residuals <- as.numeric(null_model$res)
  X1 <- as.matrix(null_model$X1)
  s2 <- as.numeric(null_model$s2)

  if (nrow(Z) != length(residuals) || nrow(X1) != nrow(Z)) {
    stop("Z and the null-model residuals have incompatible dimensions.",
         call. = FALSE)
  }
  if (length(s2) != 1L || !is.finite(s2) || s2 <= 0) {
    stop("The null model contains an invalid residual variance.",
         call. = FALSE)
  }

  K <- .modified_weighted_gaussian_kernel(Z, ell, weights)

  # This follows the modified SKAT.linear.Other calculation used for the
  # manuscript simulations.
  Q <- crossprod(residuals, K %*% residuals) / (2 * s2)
  Q_resampling <- NULL
  number_resampling <- as.integer(null_model$n.Resampling)
  if (length(number_resampling) == 1L && !is.na(number_resampling) &&
      number_resampling > 0L) {
    residuals_out <- as.matrix(null_model$res.out)
    Q_resampling <- vapply(seq_len(number_resampling), function(i) {
      as.numeric(
        crossprod(residuals_out[, i], K %*% residuals_out[, i]) /
          (2 * s2)
      )
    }, numeric(1))
  }

  XtX_inverse <- solve(crossprod(X1))
  W <- K - X1 %*% XtX_inverse %*% crossprod(X1, K)
  pvalue_matrix <- W
  if (method == "davies") {
    pvalue_matrix <- W -
      (W %*% X1) %*% XtX_inverse %*% t(X1)
  }

  pvalue_result <- .skat_mixture_pvalue(
    Q = Q,
    W = pvalue_matrix,
    Q_resampling = Q_resampling,
    method = method
  )

  fit <- list(
    p.value = as.numeric(pvalue_result$p.value),
    p.value.resampling = pvalue_result$p.value.resampling,
    Test.Type = method,
    Q = Q,
    param = pvalue_result$param,
    pval.zero.msg = pvalue_result$pval.zero.msg,
    ell = ell
  )
  fit$param$n.marker <- ncol(Z)
  fit$param$n.marker.test <- ncol(Z)
  class(fit) <- "SKAT_OUT"
  fit
}
