# GausSKAT

`GausSKAT` implements a data-adaptive rare-variant association procedure within the sequence kernel association test (SKAT) framework for continuous traits.

## Method overview

Under the default weighting scheme, let $\sqrt{w_k}=\mathrm{Beta}(\mathrm{MAF}_k;1,25)$ denote the weight assigned to variant $k$. The weighted squared genetic distance is

```math
D_{ii'}=\sum_k w_k(z_{ik}-z_{i'k})^2.
```

The weighted Gaussian kernel is

```math
K_{\ell,ii'}=\exp\left(-\frac{D_{ii'}}{\ell}\right).
```

The implementation then:

1. sets $\ell_{\min}$ to the median strictly positive weighted distance;
2. obtains the practical $\ell_{\max,\mathrm{approx}}$ from the first-order,
   covariate-projected stabilization approximation with tolerance
   $\varepsilon$;
3. constructs $J$ logarithmically spaced values between the two endpoints, with $J=5$ by default;
4. evaluates each kernel with SKAT; and
5. combines the component p-values using ACAT, with equal weights by default.

The returned object also reports the actual projected-kernel change at
$\ell_{\max,\mathrm{approx}}$, allowing the approximation to be checked directly.

## Installation

Install `GausSKAT` from GitHub:

```r
install.packages("remotes")
remotes::install_github("puneetkalra-res/GausSKAT")
```

Required dependencies, including `SKAT`, are installed automatically.

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
  number_snps = 200,
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
  seed = 2569,
  parallel = TRUE,
  n_cores = 5,
  progress = TRUE
)

simulation$summary$p_values
simulation$summary$grid
simulation$summary$endpoints
simulation$results
```

`simulation$results` contains replication-level p-values, endpoints, grid points, and endpoint diagnostics. See [`scripts/run_power_simulation.R`](scripts/run_power_simulation.R) for a command-line example.

## Reproducibility and scope

- The package currently supports continuous-trait SKAT null models.
- Genotypes must be complete dosages in `[0, 2]`; impute missing values before calling `GausSKAT()`.
- Pairwise distance and kernel matrices require $O(n^2)$ memory.
- Simulation outputs record the seed, RNG settings, and `sessionInfo()`; component-level parallelism does not alter the seeded simulation sequence.
- Dallas Heart Study participant-level data are not included in this repository.

## Relationship to SKAT

This repository adapts the continuous-trait SKAT score calculation for the weighted Gaussian kernel. It retains only the required modified pathway and uses the official [`SKAT`](https://cran.r-project.org/package=SKAT) package for null-model fitting and mixture-distribution calculations.

## Citation

If you use GausSKAT, please cite:

Kalra P, Sharma SK. GausSKAT: A Data-Adaptive Weighted Gaussian Kernel
Procedure for Rare-Variant Association Testing Within the SKAT Framework.
*Statistics in Medicine*. 2026;45(25–27):e70824.
https://doi.org/10.1002/sim.70824

Please also cite the original SKAT and ACAT publications.
