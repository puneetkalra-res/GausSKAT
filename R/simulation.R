# The stochastic call order and sampling scheme follow revision3W_simulation.R,
# the final manuscript simulation script. The SNP, haplotype, and
# causal-variant framework originates from the SKAT power code (GPL >= 2).

.validate_probability <- function(x, name, upper_inclusive = TRUE) {
  valid_upper <- if (upper_inclusive) x <= 1 else x < 1
  if (length(x) != 1L || !is.finite(x) || x < 0 || !valid_upper) {
    stop(sprintf("%s must be a valid probability.", name), call. = FALSE)
  }
  x
}

#' Define a GausSKAT simulation configuration
#'
#' @param sample_size Number of diploid individuals.
#' @param number_snps Exact number of consecutive SNPs sampled in each
#'   replication.
#' @param causal_percent Percentage of variants below `causal_maf_cutoff` that
#'   are causal.
#' @param causal_maf_cutoff Upper MAF cutoff for candidate causal variants.
#' @param negative_percent Percentage of causal effects assigned a negative
#'   sign.
#' @param effect_size_constant Constant `c` in `abs(log10(MAF)) * c`.
#' @param alpha Significance threshold used to estimate power.
#' @param weights_beta Beta-weight parameters.
#' @param epsilon GausSKAT upper-endpoint tolerance.
#' @param number_grid_points Number of GausSKAT grid points.
#' @param covariate_effects Effects of the continuous and binary covariates.
#' @param residual_sd Standard deviation of the normally distributed error.
#' @param haplotype_sampling Haplotype sampling scheme. `"manuscript"`
#'   reproduces the scheme used in the research scripts; `"full_pool"` samples
#'   from all available reference haplotypes.
#' @return A validated configuration list.
#' @export
gausskat_simulation_config <- function(
    sample_size = 1000L,
    number_snps = 200L,
    causal_percent = 50,
    causal_maf_cutoff = 0.01,
    negative_percent = 0,
    effect_size_constant = 0.15,
    alpha = 1e-6,
    weights_beta = c(1, 25),
    epsilon = 0.05,
    number_grid_points = 5L,
    covariate_effects = c(0.5, 0.5),
    residual_sd = 1,
    haplotype_sampling = c("manuscript", "full_pool")) {
  sample_size <- as.integer(sample_size)
  number_snps <- as.integer(number_snps)
  number_grid_points <- as.integer(number_grid_points)
  haplotype_sampling <- match.arg(haplotype_sampling)

  if (length(sample_size) != 1L || is.na(sample_size) || sample_size < 2L) {
    stop("sample_size must be an integer of at least two.", call. = FALSE)
  }
  if (length(number_snps) != 1L || is.na(number_snps) ||
      number_snps < 1L) {
    stop("number_snps must be a positive integer.", call. = FALSE)
  }
  if (length(causal_percent) != 1L || !is.finite(causal_percent) ||
      causal_percent <= 0 || causal_percent > 100) {
    stop("causal_percent must lie in (0, 100].", call. = FALSE)
  }
  .validate_probability(causal_maf_cutoff, "causal_maf_cutoff",
                        upper_inclusive = FALSE)
  if (causal_maf_cutoff > 0.5) {
    stop("causal_maf_cutoff cannot exceed 0.5.", call. = FALSE)
  }
  if (length(negative_percent) != 1L || !is.finite(negative_percent) ||
      negative_percent < 0 || negative_percent > 100) {
    stop("negative_percent must lie in [0, 100].", call. = FALSE)
  }
  if (length(effect_size_constant) != 1L ||
      !is.finite(effect_size_constant) || effect_size_constant < 0) {
    stop("effect_size_constant must be finite and non-negative.",
         call. = FALSE)
  }
  .validate_probability(alpha, "alpha", upper_inclusive = FALSE)
  if (alpha <= 0) {
    stop("alpha must be positive.", call. = FALSE)
  }
  if (length(weights_beta) != 2L || any(!is.finite(weights_beta)) ||
      any(weights_beta <= 0)) {
    stop("weights_beta must contain two finite positive values.",
         call. = FALSE)
  }
  if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0 ||
      epsilon >= sqrt(2)) {
    stop("epsilon must be positive and smaller than sqrt(2).",
         call. = FALSE)
  }
  if (length(number_grid_points) != 1L || is.na(number_grid_points) ||
      number_grid_points < 2L) {
    stop("number_grid_points must be at least two.", call. = FALSE)
  }
  if (length(covariate_effects) != 2L ||
      any(!is.finite(covariate_effects))) {
    stop("covariate_effects must contain two finite values.",
         call. = FALSE)
  }
  if (length(residual_sd) != 1L || !is.finite(residual_sd) ||
      residual_sd <= 0) {
    stop("residual_sd must be finite and positive.", call. = FALSE)
  }

  structure(
    list(
      sample_size = sample_size,
      number_snps = number_snps,
      causal_percent = causal_percent,
      causal_maf_cutoff = causal_maf_cutoff,
      negative_percent = negative_percent,
      effect_size_constant = effect_size_constant,
      alpha = alpha,
      weights_beta = weights_beta,
      epsilon = epsilon,
      number_grid_points = number_grid_points,
      covariate_effects = covariate_effects,
      residual_sd = residual_sd,
      haplotype_sampling = haplotype_sampling
    ),
    class = "GausSKAT_simulation_config"
  )
}

.load_skat_haplotypes <- function() {
  data_environment <- new.env(parent = emptyenv())
  utils::data("SKAT.haplotypes", package = "SKAT",
              envir = data_environment)
  if (!exists("SKAT.haplotypes", envir = data_environment,
              inherits = FALSE)) {
    stop("The SKAT.haplotypes dataset could not be loaded.", call. = FALSE)
  }
  get("SKAT.haplotypes", envir = data_environment, inherits = FALSE)
}

.prepare_haplotypes <- function(haplotypes) {
  haplotypes <- .validate_numeric_matrix(
    haplotypes,
    "haplotypes",
    minimum_rows = 2L,
    minimum_columns = 2L
  )
  if (any(haplotypes < 0 | haplotypes > 1)) {
    stop("haplotypes must contain allele indicators between 0 and 1.",
         call. = FALSE)
  }
  maf <- colMeans(haplotypes)
  # revision3W_simulation.R removes only reference-monomorphic variants with
  # allele frequency zero.
  keep <- maf > 0
  if (sum(keep) < 2L) {
    stop("Fewer than two nonzero-frequency variants remain.", call. = FALSE)
  }
  list(haplotypes = haplotypes[, keep, drop = FALSE], maf = maf[keep])
}

.sample_snps <- function(number_variants, number_snps) {
  if (number_snps >= number_variants) {
    return(seq_len(number_variants))
  }

  # revision3W_simulation.R sets SNP.dist <- seq_len(number_variants) before
  # calling Get_RandomRegion(). With consecutive integer positions, this draw
  # returns exactly number_snps consecutive SNPs (apart from probability-zero
  # integer boundary draws) while preserving the manuscript RNG sequence.
  snp_start <- stats::runif(1) *
    ((number_variants - 1) - number_snps) + 1
  snp_end <- snp_start + number_snps
  selected <- which(seq_len(number_variants) >= snp_start &
                      seq_len(number_variants) <= snp_end)
  if (length(selected) != number_snps) {
    stop("The SNP sampler did not return the requested number of SNPs.",
         call. = FALSE)
  }
  selected
}

.simulate_one_replication <- function(configuration, haplotypes, haplotype_maf,
                                      comparison_kernels, method, parallel,
                                      n_cores, worker_cluster) {
  n <- configuration$sample_size
  X1 <- stats::rnorm(n)
  X2 <- stats::rbinom(n, size = 1L, prob = 0.5)
  snps <- .sample_snps(ncol(haplotypes), configuration$number_snps)

  if (configuration$haplotype_sampling == "manuscript") {
    if (n > nrow(haplotypes)) {
      stop("sample_size exceeds the reference haplotypes under the manuscript sampling scheme.")
    }
    sampling_pool <- n
    replace_haplotypes <- n < 5000L
  } else {
    sampling_pool <- nrow(haplotypes)
    replace_haplotypes <- n > sampling_pool
  }
  # Retain the sample() expressions used by revision3W_simulation.R so the
  # seeded random-number sequence is preserved exactly.
  haplotype_1 <- sample(
    seq_len(sampling_pool),
    n,
    replace = replace_haplotypes
  )
  haplotype_2 <- sample(
    seq_len(sampling_pool),
    n,
    replace = replace_haplotypes
  )
  Z <- haplotypes[haplotype_1, snps, drop = FALSE] +
    haplotypes[haplotype_2, snps, drop = FALSE]
  region_maf <- haplotype_maf[snps]

  eligible <- which(region_maf < configuration$causal_maf_cutoff)
  if (length(eligible) == 0L) {
    stop("The sampled region contains no variants below the causal MAF cutoff.")
  }
  number_causal <- max(
    1L,
    round(configuration$causal_percent / 100 * length(eligible))
  )
  causal <- sort(sample(eligible, number_causal, replace = FALSE))
  effects <- abs(log10(region_maf[causal])) *
    configuration$effect_size_constant

  number_negative <- floor(
    length(causal) * configuration$negative_percent / 100
  )
  if (number_negative > 0L) {
    negative <- sample(seq_len(length(causal)), number_negative,
                       replace = FALSE)
    effects[negative] <- -effects[negative]
  }

  genetic_effect <- drop(Z[, causal, drop = FALSE] %*% effects)
  Y <- configuration$covariate_effects[1L] * X1 +
    configuration$covariate_effects[2L] * X2 + genetic_effect +
    stats::rnorm(n, mean = 0, sd = configuration$residual_sd)
  analysis_data <- data.frame(Y = Y, X1 = X1, X2 = X2)
  null_model <- SKAT::SKAT_Null_Model(
    Y ~ X1 + X2,
    data = analysis_data,
    out_type = "C"
  )
  X <- stats::model.matrix(~ X1 + X2, data = analysis_data)
  # Match Get_MAF_exp() in revision3W_simulation.R. The simulated reference
  # panel consists of rare variants, so these are already minor-allele
  # frequencies in the intended simulation settings.
  maf <- colMeans(Z) / 2
  weights <- .beta_weights(maf, configuration$weights_beta)

  gausskat_fit <- .gausskat_impl(
    Z = Z,
    null_model = null_model,
    X = X,
    weights = weights,
    weights_beta = configuration$weights_beta,
    epsilon = configuration$epsilon,
    number_grid_points = configuration$number_grid_points,
    method = method,
    acat_weights = NULL,
    keep_component_fits = FALSE,
    warn_on_endpoint = FALSE,
    parallel = parallel,
    n_cores = n_cores,
    worker_cluster = worker_cluster,
    call = NULL
  )

  comparison_p_values <- vapply(comparison_kernels, function(kernel) {
    fit <- SKAT::SKAT(
      Z = Z,
      obj = null_model,
      kernel = kernel,
      method = method,
      weights = weights,
      max_maf = 1
    )
    as.numeric(fit$p.value)
  }, numeric(1))

  list(
    p_gausskat = gausskat_fit$p.value,
    component_p_values = gausskat_fit$component.p.values,
    comparison_p_values = comparison_p_values,
    ell_grid = gausskat_fit$ell.grid,
    ell_min = gausskat_fit$ell.min,
    ell_max_approx = gausskat_fit$ell.max.approx,
    delta_at_ell_max = gausskat_fit$delta.at.ell.max,
    alignment_at_ell_max = gausskat_fit$alignment.at.ell.max,
    grid_degenerate = gausskat_fit$grid.degenerate,
    number_causal = length(causal)
  )
}

.simulation_summaries <- function(results, configuration,
                                  comparison_kernels) {
  method_names <- c(
    "GausSKAT",
    paste0("grid_", seq_len(configuration$number_grid_points)),
    comparison_kernels
  )
  p_columns <- c(
    "p_gausskat",
    paste0("p_grid_", seq_len(configuration$number_grid_points)),
    paste0("p_", make.names(comparison_kernels))
  )

  p_value_summary <- do.call(rbind, lapply(seq_along(p_columns), function(i) {
    values <- results[[p_columns[i]]]
    number_available <- sum(is.finite(values))
    power <- mean(values < configuration$alpha, na.rm = TRUE)
    data.frame(
      method = method_names[i],
      number_available = number_available,
      mean_p_value = mean(values, na.rm = TRUE),
      median_p_value = stats::median(values, na.rm = TRUE),
      power = power,
      power_mc_se = if (number_available > 0L && is.finite(power)) {
        sqrt(power * (1 - power) / number_available)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }))

  grid_summary <- do.call(rbind, lapply(
    seq_len(configuration$number_grid_points),
    function(j) {
      ell_values <- results[[paste0("ell_grid_", j)]]
      p_values <- results[[paste0("p_grid_", j)]]
      data.frame(
        grid_position = j,
        mean_ell = mean(ell_values, na.rm = TRUE),
        median_ell = stats::median(ell_values, na.rm = TRUE),
        mean_p_value = mean(p_values, na.rm = TRUE),
        median_p_value = stats::median(p_values, na.rm = TRUE),
        power = mean(p_values < configuration$alpha, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  ))

  grid_power <- p_value_summary[
    p_value_summary$method %in%
      paste0("grid_", seq_len(configuration$number_grid_points)),
    c("method", "power"),
    drop = FALSE
  ]
  finite_grid_power <- is.finite(grid_power$power)
  best_grid <- if (any(finite_grid_power)) {
    maximum <- max(grid_power$power[finite_grid_power])
    grid_power$method[finite_grid_power & grid_power$power == maximum]
  } else {
    character(0)
  }

  endpoint_summary <- data.frame(
    quantity = c(
      "ell_min",
      "ell_max_approx",
      "delta_at_ell_max",
      "alignment_at_ell_max"
    ),
    mean = c(
      mean(results$ell_min, na.rm = TRUE),
      mean(results$ell_max_approx, na.rm = TRUE),
      mean(results$delta_at_ell_max, na.rm = TRUE),
      mean(results$alignment_at_ell_max, na.rm = TRUE)
    ),
    median = c(
      stats::median(results$ell_min, na.rm = TRUE),
      stats::median(results$ell_max_approx, na.rm = TRUE),
      stats::median(results$delta_at_ell_max, na.rm = TRUE),
      stats::median(results$alignment_at_ell_max, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )

  list(
    p_values = p_value_summary,
    grid = grid_summary,
    endpoints = endpoint_summary,
    highest_power_grid_positions = best_grid,
    average_number_causal = mean(results$number_causal, na.rm = TRUE),
    failure_count = sum(!results$success),
    degenerate_grid_count = sum(results$grid_degenerate, na.rm = TRUE)
  )
}

#' Simulate GausSKAT power
#'
#' Runs a seed-controlled continuous-trait simulation based on the haplotypes
#' distributed with `SKAT`. Every replication records the adaptive endpoints,
#' complete grid, component p-values, aggregate p-value, and comparator
#' p-values, allowing the manuscript tables to be reconstructed from the saved
#' object.
#'
#' @param configuration A configuration returned by
#'   [gausskat_simulation_config()].
#' @param n_replications Number of simulation replications.
#' @param seed Integer random seed.
#' @param haplotypes Optional haplotype matrix. The `SKAT.haplotypes` dataset is
#'   used by default.
#' @param comparison_kernels SKAT kernels evaluated as comparators.
#' @param method SKAT p-value method.
#' @param progress Whether to display a live progress bar with percentage,
#'   elapsed time, and estimated completion time.
#' @param fail_fast Whether the first failed replication should stop the run.
#' @param parallel Whether to evaluate the component Gaussian-kernel tests in
#'   parallel within each replication. Replications themselves remain
#'   sequential so the seeded simulation sequence is unchanged.
#' @param n_cores Number of component-test workers when `parallel = TRUE`.
#'   When `NULL`, at most one fewer than the detected logical cores is used,
#'   capped by the number of grid points.
#' @return An object of class `GausSKAT_simulation`.
#' @export
simulate_gausskat_power <- function(
    configuration = gausskat_simulation_config(),
    n_replications = 100L,
    seed = 7761L,
    haplotypes = NULL,
    comparison_kernels = c("IBS.weighted", "linear.weighted"),
    method = "davies",
    progress = interactive(),
    fail_fast = TRUE,
    parallel = FALSE,
    n_cores = NULL) {
  if (!inherits(configuration, "GausSKAT_simulation_config")) {
    stop("configuration must be returned by gausskat_simulation_config().",
         call. = FALSE)
  }
  n_replications <- as.integer(n_replications)
  seed <- as.integer(seed)
  if (length(n_replications) != 1L || is.na(n_replications) ||
      n_replications < 1L) {
    stop("n_replications must be a positive integer.", call. = FALSE)
  }
  if (length(seed) != 1L || is.na(seed) || seed < 0L) {
    stop("seed must be a non-negative integer.", call. = FALSE)
  }
  if (length(comparison_kernels) < 1L ||
      any(!nzchar(comparison_kernels))) {
    stop("At least one comparison kernel must be supplied.", call. = FALSE)
  }
  if (anyDuplicated(comparison_kernels)) {
    stop("comparison_kernels must not contain duplicates.", call. = FALSE)
  }

  if (is.null(haplotypes)) {
    haplotypes <- .load_skat_haplotypes()$Haplotype
  }
  prepared <- .prepare_haplotypes(haplotypes)
  if (configuration$number_snps > ncol(prepared$haplotypes)) {
    stop("number_snps exceeds the available polymorphic variants.",
         call. = FALSE)
  }

  old_rng_kind <- RNGkind()
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv,
                            inherits = FALSE)
  if (old_seed_exists) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    do.call(RNGkind, as.list(old_rng_kind))
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  parallel_plan <- .resolve_component_parallelism(
    parallel = parallel,
    n_cores = n_cores,
    number_tasks = configuration$number_grid_points
  )
  worker_cluster <- NULL
  if (parallel_plan$backend == "PSOCK") {
    worker_cluster <- .make_psock_cluster(parallel_plan$workers)
    on.exit(parallel::stopCluster(worker_cluster), add = TRUE)
  }

  progress_bar <- NULL
  if (isTRUE(progress)) {
    progress_bar <- progress::progress_bar$new(
      format = paste0(
        "[:bar] :percent  Elapsed (h:m:s): :elapsedfull ",
        "ETA (h:m:s): :eta"
      ),
      total = n_replications,
      clear = FALSE,
      show_after = 0,
      force = TRUE
    )
    on.exit(progress_bar$terminate(), add = TRUE)
  }

  # Establish the simulation stream only after package, backend, worker, and
  # progress-bar setup. This makes the generated replications independent of
  # whether those setup operations happen to touch the RNG in a given R or
  # operating-system version. The reference manuscript script should use the
  # same boundary when exact replication-level comparison is required.
  set.seed(seed)
  manuscript_rng_kind <- RNGkind()

  runs <- vector("list", n_replications)
  for (replication in seq_len(n_replications)) {
    runs[[replication]] <- tryCatch(
      {
        value <- .simulate_one_replication(
          configuration = configuration,
          haplotypes = prepared$haplotypes,
          haplotype_maf = prepared$maf,
          comparison_kernels = comparison_kernels,
          method = method,
          parallel = parallel_plan$workers > 1L,
          n_cores = parallel_plan$workers,
          worker_cluster = worker_cluster
        )
        value$success <- TRUE
        value$error <- NA_character_
        value
      },
      error = function(error) {
        if (isTRUE(fail_fast)) {
          stop(error)
        }
        list(success = FALSE, error = conditionMessage(error))
      }
    )

    if (!is.null(progress_bar)) {
      progress_bar$tick()
    }
  }

  J <- configuration$number_grid_points
  blank_numeric <- rep(NA_real_, n_replications)
  results <- data.frame(
    replication = seq_len(n_replications),
    success = vapply(runs, function(x) isTRUE(x$success), logical(1)),
    error = vapply(runs, function(x) x$error, character(1)),
    p_gausskat = blank_numeric,
    ell_min = blank_numeric,
    ell_max_approx = blank_numeric,
    delta_at_ell_max = blank_numeric,
    alignment_at_ell_max = blank_numeric,
    grid_degenerate = rep(NA, n_replications),
    number_causal = blank_numeric,
    stringsAsFactors = FALSE
  )
  for (j in seq_len(J)) {
    results[[paste0("ell_grid_", j)]] <- blank_numeric
    results[[paste0("p_grid_", j)]] <- blank_numeric
  }
  for (kernel in comparison_kernels) {
    results[[paste0("p_", make.names(kernel))]] <- blank_numeric
  }

  for (i in which(results$success)) {
    run <- runs[[i]]
    results$p_gausskat[i] <- run$p_gausskat
    results$ell_min[i] <- run$ell_min
    results$ell_max_approx[i] <- run$ell_max_approx
    results$delta_at_ell_max[i] <- run$delta_at_ell_max
    results$alignment_at_ell_max[i] <- run$alignment_at_ell_max
    results$grid_degenerate[i] <- run$grid_degenerate
    results$number_causal[i] <- run$number_causal
    for (j in seq_len(J)) {
      results[[paste0("ell_grid_", j)]][i] <- run$ell_grid[j]
      results[[paste0("p_grid_", j)]][i] <- run$component_p_values[j]
    }
    for (kernel in comparison_kernels) {
      results[[paste0("p_", make.names(kernel))]][i] <-
        run$comparison_p_values[[kernel]]
    }
  }

  output <- list(
    configuration = configuration,
    n_replications = n_replications,
    seed = seed,
    comparison_kernels = comparison_kernels,
    method = method,
    parallel = parallel_plan$workers > 1L,
    n_cores = parallel_plan$workers,
    parallel_backend = parallel_plan$backend,
    rng_kind = manuscript_rng_kind,
    results = results,
    summary = .simulation_summaries(
      results,
      configuration,
      comparison_kernels
    ),
    session_info = utils::sessionInfo()
  )
  class(output) <- "GausSKAT_simulation"
  output
}

#' @export
print.GausSKAT_simulation <- function(x, ...) {
  cat("GausSKAT simulation\n")
  cat("  Replications: ", x$n_replications, "\n", sep = "")
  cat("  Seed:         ", x$seed, "\n", sep = "")
  cat("  Failures:     ", x$summary$failure_count, "\n", sep = "")
  cat("  Components:   ", x$parallel_backend,
      " (", x$n_cores, " worker", if (x$n_cores == 1L) "" else "s", ")\n",
      sep = "")
  print(x$summary$p_values, row.names = FALSE)
  cat("  Average number of causal variants: ",
      format(x$summary$average_number_causal, digits = 7L), "\n",
      sep = "")
  invisible(x)
}
