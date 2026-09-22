test_that("the simulation interface returns auditable replicate-level output", {
  skip_on_cran()
  skip_if_not_installed("SKAT")
  configuration <- gausskat_simulation_config(
    sample_size = 40L,
    region_length = 12L,
    causal_percent = 25,
    causal_maf_cutoff = 0.5,
    negative_percent = 20,
    effect_size_constant = 0.05,
    alpha = 0.05,
    number_grid_points = 3L,
    haplotype_sampling = "full_pool"
  )
  set.seed(98)
  haplotypes <- matrix(
    stats::rbinom(200 * 30, size = 1, prob = 0.2),
    nrow = 200
  )
  result <- simulate_gausskat_power(
    configuration,
    n_replications = 1L,
    seed = 99L,
    haplotypes = haplotypes,
    progress = FALSE,
    fail_fast = TRUE
  )

  expect_s3_class(result, "GausSKAT_simulation")
  expect_equal(result$seed, 99L)
  expect_equal(nrow(result$results), 1L)
  expect_true(result$results$success)
  expect_true(all(c("p_gausskat", "p_grid_1", "ell_grid_1") %in%
                    names(result$results)))
})
