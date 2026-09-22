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
