# GausSKAT 0.1.0

- Initial implementation of GausSKAT for continuous traits.
- Added the positive-distance median lower endpoint.
- Added the first-order, covariate-projected upper-endpoint approximation.
- Added logarithmically spaced kernel grids and equal-weight ACAT aggregation.
- Added configurable, seed-controlled power simulations.
- Added optional cross-platform parallel evaluation of Gaussian-kernel grid
  components, using forked workers on macOS/Linux and PSOCK workers on
  Windows while preserving the simulation RNG sequence.
- Integrated the focused weighted-Gaussian extension of the continuous-trait
  SKAT score calculation used in the manuscript simulations.
- Added a live simulation progress bar reporting percentage, elapsed time,
  and estimated completion time.
- Defined simulation size by the exact number of consecutive SNPs sampled in
  each replication, matching `SNP_Max` in `revision3W_simulation.R`.
- Matched the final simulation script's RNG behavior by retaining the active
  RNG kinds when applying the requested seed.
- Reported the average number of causal variants as an immediate replication-
  stream diagnostic against `revision3W_simulation.R`.
