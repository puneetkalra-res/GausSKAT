.validate_numeric_matrix <- function(x, name, minimum_rows = 1L,
                                     minimum_columns = 1L) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"

  if (nrow(x) < minimum_rows || ncol(x) < minimum_columns) {
    stop(
      sprintf(
        "%s must have at least %d row(s) and %d column(s).",
        name,
        minimum_rows,
        minimum_columns
      ),
      call. = FALSE
    )
  }

  if (any(!is.finite(x))) {
    stop(sprintf("%s contains missing or non-finite values.", name),
         call. = FALSE)
  }

  x
}

.validate_genotypes <- function(Z) {
  Z <- .validate_numeric_matrix(Z, "Z", minimum_rows = 2L)

  if (any(Z < 0 | Z > 2)) {
    stop("Z must contain complete genotype dosages between 0 and 2.",
         call. = FALSE)
  }

  Z
}

.validate_distance_matrix <- function(D) {
  D <- .validate_numeric_matrix(D, "D", minimum_rows = 2L,
                                minimum_columns = 2L)

  if (nrow(D) != ncol(D)) {
    stop("D must be a square matrix.", call. = FALSE)
  }

  scale_D <- max(1, max(abs(D)))
  symmetry_error <- max(abs(D - t(D)))

  if (symmetry_error > 1e-10 * scale_D) {
    stop("D is not symmetric to numerical tolerance.", call. = FALSE)
  }

  if (min(D) < -1e-10 * scale_D) {
    stop("D contains a materially negative value.", call. = FALSE)
  }

  D <- (D + t(D)) / 2
  D[D < 0] <- 0
  diag(D) <- 0
  D
}

.make_covariate_projection <- function(X, intercept_tolerance = 1e-10,
                                       qr_tolerance = 1e-10) {
  X <- .validate_numeric_matrix(X, "X", minimum_rows = 2L)
  qr_X <- qr(X, tol = qr_tolerance, LAPACK = FALSE)

  if (qr_X$rank < 1L) {
    stop("X has zero numerical rank.", call. = FALSE)
  }

  Q <- qr.Q(qr_X, complete = FALSE)[, seq_len(qr_X$rank), drop = FALSE]
  one <- rep(1, nrow(X))
  one_residual <- one - Q %*% drop(crossprod(Q, one))
  intercept_error <- sqrt(mean(one_residual^2))

  if (!is.finite(intercept_error) ||
      intercept_error > intercept_tolerance) {
    stop(
      paste0(
        "The column space of X must contain an intercept. RMS residual for ",
        "the intercept vector is ", signif(intercept_error, 6), "."
      ),
      call. = FALSE
    )
  }

  list(
    Q = Q,
    rank = qr_X$rank,
    number_columns = ncol(X),
    rank_deficient = qr_X$rank < ncol(X),
    intercept_error = intercept_error
  )
}

.project_symmetric_matrix <- function(M, projection) {
  M <- .validate_numeric_matrix(M, "M", minimum_rows = 2L,
                                minimum_columns = 2L)

  if (nrow(M) != ncol(M)) {
    stop("M must be a square matrix.", call. = FALSE)
  }

  Q <- projection$Q
  if (!is.matrix(Q) || nrow(Q) != nrow(M)) {
    stop("M and the covariate projection have incompatible dimensions.",
         call. = FALSE)
  }

  M_Q <- M %*% Q
  Q_M_Q <- crossprod(Q, M_Q)
  S_M_S <- M - tcrossprod(Q, M_Q) - tcrossprod(M_Q, Q) +
    tcrossprod(Q %*% Q_M_Q, Q)

  (S_M_S + t(S_M_S)) / 2
}

.included_rows <- function(Z, null_model) {
  if (is.null(null_model$id_include)) {
    return(seq_len(nrow(Z)))
  }

  id <- as.integer(null_model$id_include)
  if (length(id) < 2L || anyNA(id) || any(id < 1L) || any(id > nrow(Z))) {
    stop("The null model contains invalid included-row indices.",
         call. = FALSE)
  }
  id
}

.resolve_design_matrix <- function(X, null_model, included_rows, n_total) {
  if (is.null(X)) {
    X <- null_model$X1
  } else if (nrow(as.matrix(X)) == n_total) {
    X <- as.matrix(X)[included_rows, , drop = FALSE]
  }

  X <- .validate_numeric_matrix(X, "X", minimum_rows = 2L)
  if (nrow(X) != length(included_rows)) {
    stop("X must correspond to the observations retained by the null model.",
         call. = FALSE)
  }
  X
}

.minor_allele_frequencies <- function(Z) {
  frequency <- colMeans(Z) / 2
  pmin(frequency, 1 - frequency)
}

.beta_weights <- function(maf, weights_beta) {
  if (length(weights_beta) != 2L || any(!is.finite(weights_beta)) ||
      any(weights_beta <= 0)) {
    stop("weights_beta must contain two finite positive values.",
         call. = FALSE)
  }

  weights <- numeric(length(maf))
  polymorphic <- maf > 0 & maf < 1
  weights[polymorphic] <- stats::dbeta(
    maf[polymorphic],
    shape1 = weights_beta[1L],
    shape2 = weights_beta[2L]
  )

  if (any(!is.finite(weights))) {
    stop("The selected beta-weight parameters produced non-finite weights.",
         call. = FALSE)
  }
  weights
}
