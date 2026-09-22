#!/usr/bin/env Rscript

# Usage:
#   Rscript scripts/run_power_simulation.R [replications] [seed] [output_file] [cores]
# Example:
#   Rscript scripts/run_power_simulation.R 100 7761 results/screening.rds 5

arguments <- commandArgs(trailingOnly = TRUE)
n_replications <- if (length(arguments) >= 1L) as.integer(arguments[1L]) else 100L
seed <- if (length(arguments) >= 2L) as.integer(arguments[2L]) else 7761L
output_file <- if (length(arguments) >= 3L) {
  arguments[3L]
} else {
  file.path("results", paste0("gausskat_seed_", seed, ".rds"))
}
detected_cores <- parallel::detectCores(logical = TRUE)
if (is.na(detected_cores)) {
  detected_cores <- 1L
}
n_cores <- if (length(arguments) >= 4L) {
  as.integer(arguments[4L])
} else {
  min(5L, max(1L, detected_cores - 1L))
}
if (length(n_cores) != 1L || is.na(n_cores) || n_cores < 1L) {
  stop("cores must be a positive integer.")
}

library(GausSKAT)

configuration <- gausskat_simulation_config(
  sample_size = 1000L,
  number_snps = 200L,
  causal_percent = 50,
  causal_maf_cutoff = 0.01,
  negative_percent = 0,
  effect_size_constant = 0.15,
  alpha = 1e-6,
  weights_beta = c(1, 25),
  epsilon = 0.05,
  number_grid_points = 5L
)

result <- simulate_gausskat_power(
  configuration = configuration,
  n_replications = n_replications,
  seed = seed,
  parallel = n_cores > 1L,
  n_cores = n_cores,
  progress = TRUE,
  fail_fast = TRUE
)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(result, output_file)
write.csv(
  result$summary$p_values,
  sub("\\.rds$", "_p_value_summary.csv", output_file, ignore.case = TRUE),
  row.names = FALSE
)
write.csv(
  result$summary$endpoints,
  sub("\\.rds$", "_endpoint_summary.csv", output_file, ignore.case = TRUE),
  row.names = FALSE
)
write.csv(
  result$summary$grid,
  sub("\\.rds$", "_grid_summary.csv", output_file, ignore.case = TRUE),
  row.names = FALSE
)

print(result)
message("Saved simulation object to ", normalizePath(output_file))
