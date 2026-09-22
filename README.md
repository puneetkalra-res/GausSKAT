# GausSKAT

`GausSKAT` implements the data-adaptive weighted Gaussian-kernel sequence
kernel association test for continuous traits. It constructs a
phenotype-independent grid of weighted Gaussian kernels, obtains one SKAT
p-value per kernel, and combines the dependent p-values using the aggregated
Cauchy association test (ACAT).

This repository is an independent implementation of the method. It uses the
official [`SKAT`](https://cran.r-project.org/package=SKAT) package as a
dependency and passes each weighted Gaussian kernel through SKAT's existing
matrix-kernel interface. It does **not** redistribute or require a modified
copy of the SKAT source code.

## Method implemented

For supplied SKAT weights \(w_k\), the weighted squared genetic distance is

\[
D_{ii'}=\sum_k w_k^2(z_{ik}-z_{i'k})^2,
\]

and the weighted Gaussian kernel is

\[
K_{\ell,ii'}=\exp(-D_{ii'}/\ell).
\]

The implementation then:

1. sets \(\ell_{\min}\) to the median strictly positive weighted distance;
2. obtains the practical \(\ell_{\max,\mathrm{approx}}\) from the first-order,
   covariate-projected stabilization approximation with tolerance
   \(\varepsilon\);
3. constructs `J` logarithmically spaced values between the two endpoints;
4. evaluates each kernel with SKAT; and
5. combines the component p-values by equal-weight ACAT.

The returned object reports the actual projected-kernel change at
\(\ell_{\max,\mathrm{approx}}\), so the approximation is transparent and can
be checked for every analysis or simulation replication.

## Installation

Install the official dependency and then install this repository:

```r
install.packages("SKAT")
install.packages("remotes")
remotes::install_github("puneetkalra-res/GausSKAT")
```

For a local checkout:

```r
install.packages("devtools")
devtools::install(".")
```

## Basic analysis

```r
library(SKAT)
library(GausSKAT)

data("SKAT.example", package = "SKAT")
null_model <- SKAT_Null_Model(
  y.c ~ X,
  data = SKAT.example,
  out_type = "C"
)

fit <- GausSKAT(
  Z = SKAT.example$Z,
  null_model = null_model,
  epsilon = 0.05,
  number_grid_points = 5,
  parallel = TRUE,
  n_cores = 5
)

fit
fit$p.value
fit$component.p.values
fit$ell.grid
```

The grid is constructed from the genotype matrix, variant weights, and the
null-model covariate design; it does not use the phenotype values.

## Power simulation

```r
configuration <- gausskat_simulation_config(
  sample_size = 1000,
  region_length = 200,
  causal_percent = 50,
  causal_maf_cutoff = 0.01,
  negative_percent = 0,
  effect_size_constant = 0.15,
  alpha = 1e-6,
  epsilon = 0.05,
  number_grid_points = 5
)

simulation <- simulate_gausskat_power(
  configuration = configuration,
  n_replications = 100,
  seed = 7761,
  parallel = TRUE,
  n_cores = 5,
  progress = TRUE
)

simulation$summary$p_values
simulation$summary$grid
simulation$summary$endpoints
simulation$results
```

`simulation$results` retains every p-value, endpoint, grid point, and endpoint
diagnostic needed to reproduce detailed supplementary tables. See
[`scripts/run_power_simulation.R`](scripts/run_power_simulation.R) for a
command-line example.

## Reproducibility and scope

- The package currently supports continuous-trait SKAT null models.
- Genotypes must be complete dosages in `[0, 2]`; perform and document any
  imputation before calling `GausSKAT()`.
- The method allocates pairwise distance and kernel matrices and therefore has
  \(O(n^2)\) memory use.
- The first-order upper endpoint is labelled `ell.max.approx` throughout the
  code and output. Its observed doubling-change diagnostic is also returned.
- Simulation runs explicitly use `Mersenne-Twister`, `Inversion`, and
  `Rejection`, matching `revision3W_simulation.R` when run in a fresh R session.
  The seed, RNG settings, and `sessionInfo()` are recorded, and the caller's
  pre-existing random-number state is restored on exit.
- Optional component-level parallelism uses forked workers on macOS/Linux and
  a persistent PSOCK cluster on Windows. Simulation replications remain
  sequential, and only deterministic kernel tests are parallelised, so the
  seeded simulation sequence is unchanged.
- The default `haplotype_sampling = "manuscript"` reproduces the sampling
  scheme in the research scripts. The alternative `"full_pool"` scheme is
  available for new experiments and is recorded in the configuration.
- No Dallas Heart Study participant-level data, intermediate results, or local
  machine paths are included in this repository.

## Repository layout

- `R/`: method and simulation functions.
- `tests/testthat/`: numerical, interface, and smoke tests.
- `scripts/`: executable analysis and reproducibility examples.
- `vignettes/`: a concise method walkthrough.
- `.github/workflows/`: automated package checks.
- `CODE_PROVENANCE.md`: separation of the final implementation from the
  exploratory modified SKAT source.
- `PUBLISHING.md`: the short sequence for creating and pushing the repository.

## Relationship to SKAT

The SKAT authors retain ownership of the `SKAT` package. `GausSKAT` calls the
installed package and relies on its score-test and mixture-distribution
calculations. The external-kernel construction is tested against the weighted
Gaussian formula used in the research code; see
[`scripts/verify_kernel_equivalence.R`](scripts/verify_kernel_equivalence.R).

## Citation

Please cite the accompanying GausSKAT manuscript and the original SKAT and
ACAT publications. Full manuscript citation details can be added to
`CITATION.cff` when the article DOI is assigned.
