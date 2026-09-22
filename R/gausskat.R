.resolve_component_parallelism <- function(parallel, n_cores,
                                           number_tasks) {
  if (length(parallel) != 1L || is.na(parallel) || !is.logical(parallel)) {
    stop("parallel must be TRUE or FALSE.", call. = FALSE)
  }

  number_tasks <- as.integer(number_tasks)
  if (length(number_tasks) != 1L || is.na(number_tasks) ||
      number_tasks < 1L) {
    stop("number_tasks must be a positive integer.", call. = FALSE)
  }

  if (is.null(n_cores)) {
    detected <- parallel::detectCores(logical = TRUE)
    if (is.na(detected)) {
      detected <- 1L
    }
    n_cores <- max(1L, detected - 1L)
  } else {
    n_cores <- as.integer(n_cores)
    if (length(n_cores) != 1L || is.na(n_cores) || n_cores < 1L) {
      stop("n_cores must be NULL or a positive integer.", call. = FALSE)
    }
  }

  workers <- if (isTRUE(parallel)) {
    min(n_cores, number_tasks)
  } else {
    1L
  }
  backend <- if (workers == 1L) {
    "sequential"
  } else if (.Platform$OS.type == "windows") {
    "PSOCK"
  } else {
    "fork"
  }

  list(workers = workers, backend = backend)
}

.make_psock_cluster <- function(n_cores) {
  seed_exists <- exists(".Random.seed", envir = .GlobalEnv,
                        inherits = FALSE)
  if (seed_exists) {
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (seed_exists) {
      assign(".Random.seed", seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv,
                      inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  parallel::makeCluster(n_cores, type = "PSOCK", useXDR = FALSE)
}

.evaluate_component_fits <- function(ell_grid, Z, null_model, weights, method,
                                     parallel_plan,
                                     worker_cluster = NULL) {
  evaluate_one <- function(ell) {
    .modified_skat_gaussian(
      Z = Z,
      null_model = null_model,
      ell = ell,
      weights = weights,
      method = method
    )
  }

  if (parallel_plan$backend == "sequential") {
    return(lapply(ell_grid, evaluate_one))
  }

  if (parallel_plan$backend == "fork") {
    return(parallel::mclapply(
      ell_grid,
      evaluate_one,
      mc.cores = parallel_plan$workers,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    ))
  }

  cluster_created_here <- is.null(worker_cluster)
  if (cluster_created_here) {
    worker_cluster <- .make_psock_cluster(parallel_plan$workers)
    on.exit(parallel::stopCluster(worker_cluster), add = TRUE)
  }

  parallel::clusterExport(
    worker_cluster,
    varlist = c("Z", "null_model", "weights", "method"),
    envir = environment()
  )

  worker_function <- function(ell) {
    fit_function <- getFromNamespace(
      ".modified_skat_gaussian",
      "GausSKAT"
    )
    fit <- fit_function(
      Z = Z,
      null_model = null_model,
      ell = ell,
      weights = weights,
      method = method
    )
    if (length(fit$p.value) != 1L || !is.finite(fit$p.value)) {
      stop("SKAT did not return a finite component p-value.",
           call. = FALSE)
    }
    fit
  }
  # The required matrices are already present in each worker's global
  # environment. Avoid serializing the enclosing function environment, which
  # would otherwise transmit the same large objects a second time.
  environment(worker_function) <- .GlobalEnv

  unname(parallel::parLapplyLB(
    worker_cluster,
    ell_grid,
    worker_function
  ))
}

#' GausSKAT test for a continuous trait
#'
#' Constructs a phenotype-independent, data-adaptive grid of weighted Gaussian
#' kernels, evaluates each kernel using the weighted-Gaussian extension of the
#' continuous-trait SKAT score calculation, and
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
#' @param method SKAT p-value method.
#' @param acat_weights Optional ACAT component weights. Equal weights are used
#'   by default.
#' @param keep_component_fits Whether to retain the complete SKAT fit for every
#'   grid point.
#' @param warn_on_endpoint Whether to warn when the observed doubling change at
#'   the approximate upper endpoint exceeds `epsilon`.
#' @param parallel Whether to evaluate the component Gaussian-kernel tests in
#'   parallel. The default is `FALSE`.
#' @param n_cores Number of component-test workers when `parallel = TRUE`.
#'   When `NULL`, at most one fewer than the detected logical cores is used,
#'   capped by the number of grid points.
#' @return An object of class `GausSKAT` containing the aggregated p-value,
#'   component p-values, adaptive grid, and endpoint diagnostics.
#' @export
GausSKAT <- function(Z, null_model, X = NULL, weights = NULL,
                     weights_beta = c(1, 25), epsilon = 0.05,
                     number_grid_points = 5L, method = "davies",
                     acat_weights = NULL, keep_component_fits = FALSE,
                     warn_on_endpoint = TRUE, parallel = FALSE,
                     n_cores = NULL) {
  .gausskat_impl(
    Z = Z,
    null_model = null_model,
    X = X,
    weights = weights,
    weights_beta = weights_beta,
    epsilon = epsilon,
    number_grid_points = number_grid_points,
    method = method,
    acat_weights = acat_weights,
    keep_component_fits = keep_component_fits,
    warn_on_endpoint = warn_on_endpoint,
    parallel = parallel,
    n_cores = n_cores,
    worker_cluster = NULL,
    call = match.call()
  )
}

.gausskat_impl <- function(Z, null_model, X, weights, weights_beta, epsilon,
                           number_grid_points, method, acat_weights,
                           keep_component_fits, warn_on_endpoint, parallel,
                           n_cores, worker_cluster, call) {
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

  parallel_plan <- .resolve_component_parallelism(
    parallel = parallel,
    n_cores = n_cores,
    number_tasks = number_grid_points
  )

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

  component_fits <- .evaluate_component_fits(
    ell_grid = grid$ell_grid,
    Z = Z_analysis,
    null_model = null_model,
    weights = weights,
    method = method,
    parallel_plan = parallel_plan,
    worker_cluster = worker_cluster
  )
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
    parallel = parallel_plan$workers > 1L,
    n.cores = parallel_plan$workers,
    parallel.backend = parallel_plan$backend,
    computation.path = "modified weighted-Gaussian SKAT",
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
  cat("  Component backend:  ", x$parallel.backend,
      " (", x$n.cores, " worker", if (x$n.cores == 1L) "" else "s", ")\n",
      sep = "")
  invisible(x)
}
