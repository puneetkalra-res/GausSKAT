# Code provenance and separation from SKAT

The research prototype was developed by importing the SKAT source into an
RStudio project and adding a weighted Gaussian kernel with a user-supplied
hyperparameter. The three changes essential to that prototype were:

1. construction of `weighted_Z` by multiplying genotype column `k` by the
   supplied SKAT weight;
2. construction of `K = exp(-dist(weighted_Z)^2 / ell)`; and
3. propagation of `ell` through the local SKAT call stack.

This repository incorporates the continuous-trait weighted-Gaussian pathway
used in the final research code. The file `R/modified_skat_gaussian.R`
isolates the weighted-Gaussian kernel construction and the corresponding
continuous-trait calculation adapted from `SKAT.linear.Other`. The installed
`SKAT` package remains the dependency for null-model fitting and the
established Davies and Liu mixture-distribution calculations. This avoids the
slower external matrix-kernel route while retaining only the source changes
needed by GausSKAT.

The repository deliberately excludes all other changes present in the local
research copy, including alternative binary-trait links, GPU experiments,
package renaming, altered default kernels, obsolete endpoint searches, and
local compiled objects. It also excludes all participant-level Dallas Heart
Study files and intermediate simulation output.

The script `scripts/verify_kernel_equivalence.R` compares the package
kernel constructor with the prototype formula on a fixed example. Automated
tests additionally compare the modified pathway against the matrix-kernel
calculation, and cover endpoint construction, weight-scale equivariance, the
large-hyperparameter weighted-linear limit, ACAT, and the simulation output
structure.

The default simulation RNG sequence follows `revision3W_simulation.R`, the
final full-simulation script: a fresh-session Mersenne--Twister stream is seeded
once and used continuously across replications. The normal and sampling
generators are fixed explicitly as `Inversion` and `Rejection` so that prior
session settings cannot silently change the sequence.

Optional parallel execution is restricted to the deterministic Gaussian-kernel
component tests within each replication. Simulation replications retain their
original sequential order. The implementation uses forked workers on
macOS/Linux and a persistent PSOCK cluster on Windows, without advancing the
master random-number stream.
